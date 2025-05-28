#!/usr/bin/env python3
"""
Genome Downloader
Download and organize bacterial genomes from NCBI using assembly_summary.txt
"""

import os
import sys
import subprocess
import glob
import zipfile
import shutil
from datetime import datetime
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import time
import pandas as pd


class GenomeDownloader:
    """
    Download and organize bacterial genomes from NCBI using assembly_summary.txt
    """
    
    def __init__(self, work_dir=None, assembly_summary_path=None):
        """
        Initialize the genome downloader
        
        Args:
            work_dir (str): Working directory (default: current directory)
            assembly_summary_path (str): Path to assembly_summary.txt file
        """
        self.work_dir = Path(work_dir) if work_dir else Path.cwd()
        
        # Simplified assembly_summary.txt path detection
        if assembly_summary_path:
            self.assembly_summary_path = assembly_summary_path
        else:
            # Default to data directory in the project
            default_path = Path(__file__).parent.parent / "data" / "assembly_summary.txt"
            self.assembly_summary_path = str(default_path)
        
        # Check if assembly_summary.txt exists
        if not os.path.exists(self.assembly_summary_path):
            print(f"❌ Error: Cannot find {self.assembly_summary_path}")
            print("Please ensure assembly_summary.txt is available in src/microbial_primer_design/data/")
            print("You can download it with:")
            print("wget https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt -O src/microbial_primer_design/data/assembly_summary.txt")
            raise FileNotFoundError(f"assembly_summary.txt not found at {self.assembly_summary_path}")
        
        print(f"📁 Using assembly_summary.txt at: {self.assembly_summary_path}")
    
    def extract_genomes_by_level(self, genus, level="Complete Genome"):
        """
        Extract genome accessions for a specific genus from assembly_summary.txt
        
        Args:
            genus (str): Target genus name
            level (str): Assembly level filter
            
        Returns:
            str: Path to the generated list file, or None if no genomes found
        """
        genus_dir = self.work_dir / genus
        genus_dir.mkdir(exist_ok=True)
        
        list_file = genus_dir / f"{genus}_list.txt"
        found = False

        with open(self.assembly_summary_path, "r") as infile, open(list_file, "w") as outfile:
            for line in infile:
                if genus in line and level in line:
                    outfile.write(line)
                    found = True

        return str(list_file) if found else None
    
    def extract_genome_summary(self, genus, level="Complete Genome"):
        """
        Extract detailed genome summary information for a specific genus
        
        Args:
            genus (str): Target genus name
            level (str): Assembly level filter
            
        Returns:
            str: Path to the generated summary file, or None if no genomes found
        """
        genus_dir = self.work_dir / genus
        genus_dir.mkdir(exist_ok=True)
        
        summary_file = genus_dir / f"{genus}_genome_summary.csv"
        
        # Read assembly_summary.txt and filter for target genus
        try:
            # Read the header first
            with open(self.assembly_summary_path, "r") as f:
                header = f.readline().strip()
            
            # Read the full file
            df = pd.read_csv(self.assembly_summary_path, sep='\t', low_memory=False, comment='#')
            
            # Filter for target genus and assembly level
            filtered_df = df[
                (df['organism_name'].str.contains(genus, case=False, na=False)) &
                (df['assembly_level'] == level)
            ]
            
            if len(filtered_df) == 0:
                print(f"⚠️  No genomes found for {genus} with level '{level}'")
                return None
            
            # Select important columns for summary
            summary_columns = [
                'assembly_accession',
                'bioproject',
                'biosample', 
                'wgs_master',
                'refseq_category',
                'taxid',
                'species_taxid',
                'organism_name',
                'infraspecific_name',
                'isolate',
                'version_status',
                'assembly_level',
                'release_type',
                'genome_rep',
                'seq_rel_date',
                'asm_name',
                'submitter',
                'gbrs_paired_asm',
                'paired_asm_comp',
                'ftp_path',
                'excluded_from_refseq',
                'relation_to_type_material'
            ]
            
            # Keep only available columns
            available_columns = [col for col in summary_columns if col in filtered_df.columns]
            summary_df = filtered_df[available_columns].copy()
            
            # Save summary
            summary_df.to_csv(summary_file, index=False)
            
            print(f"📊 Extracted summary for {len(summary_df)} {genus} genomes")
            print(f"📁 Summary saved to: {summary_file}")
            
            # Display basic statistics
            print(f"📈 Summary statistics:")
            print(f"   Total genomes: {len(summary_df)}")
            if 'refseq_category' in summary_df.columns:
                refseq_counts = summary_df['refseq_category'].value_counts()
                for category, count in refseq_counts.items():
                    print(f"   {category}: {count}")
            
            return str(summary_file)
            
        except Exception as e:
            print(f"❌ Failed to extract genome summary: {e}")
            return None
    
    def download_single_genome(self, accession):
        """
        Download a single genome by accession number
        
        Args:
            accession (str): Genome accession number
            
        Returns:
            bool: Success status
        """
        filename = f"{accession}.zip"
        try:
            print(f"🔄 Downloading genome {accession}...")
            result = subprocess.run(
                ['datasets', 'download', 'genome', 'accession', accession, '--filename', filename],
                check=True,
                capture_output=True,
                text=True,
                cwd=self.work_dir
            )
            print(f"✅ Successfully downloaded genome {accession}")
            return True
        except subprocess.CalledProcessError as e:
            print(f"❌ Failed to download genome {accession}: {e.stderr}")
            return False
    
    def batch_download_genomes(self, listfile, batch_size=10):
        """
        Download genomes in batches from a list file
        
        Args:
            listfile (str): Path to file containing accession numbers
            batch_size (int): Number of concurrent downloads
            
        Returns:
            bool: Success status
        """
        # Read all accession numbers
        with open(listfile, 'r') as f:
            accessions = [line.strip().split('\t')[0] for line in f.readlines()]
        
        total_accessions = len(accessions)
        print(f"📊 Total genomes to download: {total_accessions}")
        
        success_count = 0
        
        # Use thread pool for parallel downloads
        with ThreadPoolExecutor(max_workers=batch_size) as executor:
            for i in range(0, total_accessions, batch_size):
                batch = accessions[i:i + batch_size]
                print(f"\n🚀 Processing batch {i//batch_size + 1}, {len(batch)} genomes")
                
                # Submit all download tasks for current batch
                futures = [executor.submit(self.download_single_genome, acc) for acc in batch]
                
                # Wait for current batch to complete
                for future in futures:
                    if future.result():
                        success_count += 1
                
                # Pause between batches to avoid server overload
                if i + batch_size < total_accessions:
                    print("⏳ Waiting 5 seconds before next batch...")
                    time.sleep(5)
        
        print(f"📈 Download completed: {success_count}/{total_accessions} successful")
        return success_count > 0
    
    def unzip_files(self, genus_dir):
        """
        Extract all zip files in the specified directory
        
        Args:
            genus_dir (Path): Directory containing zip files
            
        Returns:
            bool: Success status
        """
        print(f"🔄 [{genus_dir.name}] Extracting zip files...")
        
        zip_files = list(genus_dir.glob("*.zip"))
        if not zip_files:
            print(f"⚠️  [{genus_dir.name}] No zip files found")
            return False
        
        for zip_file in zip_files:
            try:
                # Create extraction directory (remove .zip extension)
                extract_dir = zip_file.with_suffix('')
                with zipfile.ZipFile(zip_file, 'r') as zip_ref:
                    zip_ref.extractall(extract_dir)
                print(f"✅ Extracted: {zip_file.name}")
            except Exception as e:
                print(f"❌ Failed to extract {zip_file.name}: {e}")
                return False
        
        return True
    
    def collect_fna_files(self, genus_dir):
        """
        Collect all fna files to the genus directory root
        
        Args:
            genus_dir (Path): Genus directory
            
        Returns:
            bool: Success status
        """
        print(f"🔄 [{genus_dir.name}] Collecting fna files...")
        
        # Find all fna files (may be in multiple subdirectories)
        fna_files = list(genus_dir.rglob("*.fna"))
        
        if not fna_files:
            print(f"⚠️  [{genus_dir.name}] No fna files found")
            return False
        
        # Move all fna files to genus directory root
        for fna_file in fna_files:
            if fna_file.parent != genus_dir:  # Only move files not in root
                dest_file = genus_dir / fna_file.name
                try:
                    shutil.move(str(fna_file), str(dest_file))
                    print(f"✅ Moved file: {fna_file.name}")
                except Exception as e:
                    print(f"❌ Failed to move file {fna_file.name}: {e}")
        
        return True
    
    def organize_genus_files(self, genus, is_outgroup=False):
        """
        Organize genus files into final directory structure
        
        Args:
            genus (str): Genus name
            is_outgroup (bool): Whether this is an outgroup genus
            
        Returns:
            bool: Success status
        """
        genus_dir = self.work_dir / genus
        
        if not genus_dir.exists():
            print(f"❌ Directory does not exist: {genus_dir}")
            return False
        
        print(f"🔄 [{genus}] Organizing file structure...")
        
        # Create target directory
        if is_outgroup:
            target_base_dir = genus_dir
            target_data_dir = target_base_dir / "data" / "outgroup"
            prefix = "Out"
        else:
            target_base_dir = genus_dir
            target_data_dir = target_base_dir / "data" / genus[:3].capitalize()
            prefix = genus[:3].capitalize()
        
        target_data_dir.mkdir(parents=True, exist_ok=True)
        
        # Get all fna files
        fna_files = list(genus_dir.glob("*.fna"))
        
        if not fna_files:
            print(f"⚠️  [{genus}] No fna files found for organization")
            return False
        
        # Rename and move fna files
        for fna_file in fna_files:
            original_name = fna_file.name
            new_name = f"{prefix}_{original_name}"
            target_file = target_data_dir / new_name
            
            try:
                shutil.move(str(fna_file), str(target_file))
                print(f"✅ Organized file: {original_name} -> {new_name}")
            except Exception as e:
                print(f"❌ Failed to organize file {original_name}: {e}")
        
        print(f"✅ [{genus}] File organization completed, saved to: {target_data_dir}")
        return True
    
    def process_genus_complete(self, genus, level="Complete Genome", is_outgroup=False):
        """
        Complete processing of a genus: download, extract, organize files
        
        Args:
            genus (str): Genus name
            level (str): Assembly level filter
            is_outgroup (bool): Whether this is an outgroup genus
            
        Returns:
            tuple: (genus_name, success_status)
        """
        print(f"\n🚀 Starting processing {'outgroup ' if is_outgroup else ''}genus: {genus}")
        
        # Step 1: Extract genome list and download
        if not self.download_genus(genus, level):
            return genus, False
        
        genus_dir = self.work_dir / genus
        
        # Step 2: Extract zip files
        if not self.unzip_files(genus_dir):
            return genus, False
        
        # Step 3: Collect fna files
        if not self.collect_fna_files(genus_dir):
            return genus, False
        
        # Step 4: Organize file structure
        if not self.organize_genus_files(genus, is_outgroup):
            return genus, False
        
        print(f"🎉 [{genus}] Complete processing finished!")
        return genus, True
    
    def download_genus(self, genus, level="Complete Genome"):
        """
        Download genomes for a specific genus and extract summary information
        
        Args:
            genus (str): Target genus name
            level (str): Assembly level filter
            
        Returns:
            bool: Success status
        """
        genus_dir = self.work_dir / genus
        genus_dir.mkdir(exist_ok=True)
        log_file = genus_dir / "download.log"
        not_found_log = self.work_dir / "not_found.log"

        # Extract genome list from assembly_summary.txt
        list_file = self.extract_genomes_by_level(genus, level)
        if not list_file:
            with open(not_found_log, "a") as nf:
                nf.write(f"[{datetime.now()}] NOT FOUND: {genus} (level: {level})\n")
            print(f"❌ No genomes found for {genus} with level '{level}'")
            return False
        
        # Extract detailed genome summary information
        print(f"📊 Extracting genome summary information for {genus}...")
        summary_file = self.extract_genome_summary(genus, level)
        if summary_file:
            print(f"✅ Genome summary extracted successfully")
        else:
            print(f"⚠️  Failed to extract genome summary, but continuing with download")
        
        # Download genomes
        try:
            success = self.batch_download_genomes(list_file)
            with open(log_file, "a") as f:
                status = "SUCCESS" if success else "FAIL"
                f.write(f"[{datetime.now()}] {status}: {genus}\n")
                if summary_file:
                    f.write(f"[{datetime.now()}] SUMMARY: {summary_file}\n")
            return success
        except Exception as e:
            with open(log_file, "a") as f:
                f.write(f"[{datetime.now()}] ERROR: {genus} - {e}\n")
            print(f"❌ Error downloading {genus}: {e}")
            return False
    
    def download_with_outgroup(self, target_genus, outgroup_genera, 
                              level="Complete Genome", threads=4):
        """
        Download target genus and outgroup genomes
        
        Args:
            target_genus (str): Target genus name
            outgroup_genera (list): List of outgroup genus names
            level (str): Assembly level filter
            threads (int): Number of parallel threads
            
        Returns:
            bool: Success status
        """
        print(f"⚙️ Starting parallel download and processing, threads: {threads}")
        print(f"📋 Target genera: {target_genus}")
        if outgroup_genera:
            print(f"🔗 Outgroups: {', '.join(outgroup_genera)}")

        all_genera = [target_genus] + (outgroup_genera or [])

        # Use thread pool for concurrent processing
        with ThreadPoolExecutor(max_workers=threads) as executor:
            futures = {}
            
            # Submit target genus task
            future = executor.submit(self.process_genus_complete, target_genus, level, False)
            futures[future] = (target_genus, False)
            
            # Submit outgroup tasks
            for outgroup in outgroup_genera or []:
                future = executor.submit(self.process_genus_complete, outgroup, level, True)
                futures[future] = (outgroup, True)

            # Process results
            success_count = 0
            for future in as_completed(futures):
                genus, is_outgroup = futures[future]
                try:
                    genus_result, success = future.result()
                    group_type = "outgroup" if is_outgroup else "target"
                    if success:
                        print(f"🎯 [{group_type}] {genus} processing completed successfully")
                        success_count += 1
                    else:
                        print(f"⚠️  [{group_type}] {genus} processing failed")
                except Exception as exc:
                    group_type = "outgroup" if is_outgroup else "target"
                    print(f"❌ [{group_type}] {genus} execution error: {exc}")

        print("\n🏁 All tasks completed!")
        print("📁 File organization structure:")
        
        # Show final directory structure
        target_data_dir = self.work_dir / target_genus / "data" / target_genus[:3].capitalize()
        if target_data_dir.exists():
            print(f"   {target_data_dir}/")
        
        for outgroup in outgroup_genera or []:
            outgroup_data_dir = self.work_dir / outgroup / "data" / "outgroup"
            if outgroup_data_dir.exists():
                print(f"   {outgroup_data_dir}/")

        return success_count > 0


def main():
    """
    Command line interface for GenomeDownloader
    """
    import argparse
    
    parser = argparse.ArgumentParser(description="Download bacterial genomes using assembly_summary.txt")
    parser.add_argument("genus", help="Target genus name")
    parser.add_argument("--outgroup", nargs="*", help="Outgroup genus names")
    parser.add_argument("--level", default="Complete Genome", 
                       help="Assembly level (default: 'Complete Genome')")
    parser.add_argument("--threads", type=int, default=4, 
                       help="Number of threads (default: 4)")
    parser.add_argument("--work-dir", help="Working directory")
    parser.add_argument("--assembly-summary", help="Path to assembly_summary.txt file")
    
    args = parser.parse_args()
    
    try:
        # Create downloader
        downloader = GenomeDownloader(args.work_dir, args.assembly_summary)
        
        # Download genomes
        if args.outgroup:
            success = downloader.download_with_outgroup(
                args.genus, args.outgroup, args.level, args.threads
            )
        else:
            success = downloader.download_genus(args.genus, args.level)
        
        if success:
            print("✅ Download completed successfully")
        else:
            print("❌ Download failed")
            sys.exit(1)
            
    except FileNotFoundError as e:
        print(f"❌ {e}")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Unexpected error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main() 