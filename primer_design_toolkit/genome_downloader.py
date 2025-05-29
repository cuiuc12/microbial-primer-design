"""
Genome Downloader
Download and organize bacterial genomes from NCBI using assembly_summary.txt
Enhanced with improved file organization strategy to prevent target/outgroup confusion
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
    
    File organization strategy:
    - Target genus: target_genus/data/Target_genus_prefix/
    - Outgroups: target_genus/data/outgroup/
    """
    
    def __init__(self, work_dir=None, assembly_summary_path=None):
        """
        初始化基因组下载器 / Initialize the genome downloader
        
        Args:
            work_dir (str): 工作目录 / Working directory (default: current directory)
            assembly_summary_path (str): assembly_summary.txt文件路径 / Path to assembly_summary.txt file
        """
        self.work_dir = Path(work_dir) if work_dir else Path.cwd()
        
        # 简化assembly_summary.txt路径检测 / Simplified assembly_summary.txt path detection
        if assembly_summary_path:
            self.assembly_summary_path = assembly_summary_path
        else:
            # 默认到项目的data目录 / Default to data directory in the project
            default_path = Path(__file__).parent.parent / "data" / "assembly_summary.txt"
            self.assembly_summary_path = str(default_path)
        
        # 检查assembly_summary.txt是否存在 / Check if assembly_summary.txt exists
        if not os.path.exists(self.assembly_summary_path):
            print(f"❌ Error: Cannot find {self.assembly_summary_path}")
            print("Please ensure assembly_summary.txt is available in src/microbial_primer_design/data/")
            print("You can download it with:")
            print("wget https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt -O src/microbial_primer_design/data/assembly_summary.txt")
            raise FileNotFoundError(f"assembly_summary.txt not found at {self.assembly_summary_path}")
        
        print(f"📁 Using assembly_summary.txt at: {self.assembly_summary_path}")
    
    def extract_genomes_by_level(self, genus, level="Complete Genome"):
        """
        提取指定genus的基因组accession / Extract genome accessions for a specific genus
        Args:
            genus (str): 目标属名 / Target genus name
            level (str): 组装级别过滤 / Assembly level filter
        Returns:
            str: 生成的list文件路径 / Path to the generated list file
        """
        genus_safe = genus.replace(" ", "_")  # 用于文件/目录名的安全属名 / Safe genus name for file/dir
        genus_dir = self.work_dir / genus_safe  # 属目录 / Genus directory
        genus_dir.mkdir(exist_ok=True)
        list_file = genus_dir / f"{genus_safe}_list.txt"
        found = False
        with open(self.assembly_summary_path, "r") as infile, open(list_file, "w") as outfile:
            for line in infile:
                if genus in line and level in line:  # 用原始genus做过滤 / Use original genus for filtering
                    outfile.write(line)
                    found = True
        return str(list_file) if found else None
    
    def extract_genome_summary(self, genus, level="Complete Genome"):
        """
        提取指定genus的详细summary信息 / Extract detailed genome summary information for a specific genus
        Args:
            genus (str): 目标属名 / Target genus name
            level (str): 组装级别过滤 / Assembly level filter
        Returns:
            str: 生成的summary文件路径 / Path to the generated summary file
        """
        genus_safe = genus.replace(" ", "_")  # 用于文件/目录名的安全属名 / Safe genus name for file/dir
        genus_dir = self.work_dir / genus_safe  # 属目录 / Genus directory
        genus_dir.mkdir(exist_ok=True)
        summary_file = genus_dir / f"{genus_safe}_genome_summary.csv"
        try:
            # 精确查找表头行号和列名 / Find header line number and column names
            header_line = None
            header_cols = None
            with open(self.assembly_summary_path, "r") as f:
                for i, line in enumerate(f):
                    if line.lstrip().startswith("#assembly_accession"):
                        header_line = i
                        header_cols = [col.lstrip('#').strip() for col in line.strip().split('\t')]
                        break
            if header_line is None or header_cols is None:
                print(f"❌ Cannot find header line in {self.assembly_summary_path}")
                return None
            # 用skiprows跳过表头前的注释行，并手动指定列名 / Use skiprows and names to set column names
            df = pd.read_csv(
                self.assembly_summary_path,
                sep='\t',
                skiprows=header_line + 1,
                names=header_cols,
                low_memory=False
            )
            filtered_df = df[
                (df['organism_name'].str.contains(genus, case=False, na=False)) &
                (df['assembly_level'] == level)
            ]
            if len(filtered_df) == 0:
                print(f"⚠️  No genomes found for {genus} with level '{level}'")
                return None
            summary_columns = [
                'assembly_accession', 'bioproject', 'biosample', 'wgs_master', 'refseq_category',
                'taxid', 'species_taxid', 'organism_name', 'infraspecific_name', 'isolate',
                'version_status', 'assembly_level', 'release_type', 'genome_rep', 'seq_rel_date',
                'asm_name', 'submitter', 'gbrs_paired_asm', 'paired_asm_comp', 'ftp_path',
                'excluded_from_refseq', 'relation_to_type_material'
            ]
            available_columns = [col for col in summary_columns if col in filtered_df.columns]
            summary_df = filtered_df[available_columns].copy()
            summary_df.to_csv(summary_file, index=False)
            print(f"📊 Extracted summary for {len(summary_df)} {genus} genomes")
            print(f"📁 Summary saved to: {summary_file}")
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
    
    def extract_genome_summary_to_dir(self, genus, level="Complete Genome", target_dir=None):
        """
        提取指定genus的详细summary信息到指定目录 / Extract detailed genome summary information for a specific genus to specified directory
        Args:
            genus (str): 目标属名 / Target genus name
            level (str): 组装级别过滤 / Assembly level filter
            target_dir (Path): 目标保存目录 / Target save directory
        Returns:
            str: 生成的summary文件路径 / Path to the generated summary file
        """
        genus_safe = genus.replace(" ", "_")
        
        # 确定保存目录 / Determine save directory
        if target_dir is not None:
            save_dir = target_dir
        else:
            save_dir = self.work_dir / genus_safe
        save_dir.mkdir(exist_ok=True)
        
        summary_file = save_dir / f"{genus_safe}_genome_summary.csv"
        
        try:
            # 精确查找表头行号和列名 / Find header line number and column names
            header_line = None
            header_cols = None
            with open(self.assembly_summary_path, "r") as f:
                for i, line in enumerate(f):
                    if line.lstrip().startswith("#assembly_accession"):
                        header_line = i
                        header_cols = [col.lstrip('#').strip() for col in line.strip().split('\t')]
                        break
            if header_line is None or header_cols is None:
                print(f"❌ Cannot find header line in {self.assembly_summary_path}")
                return None
            # 用skiprows跳过表头前的注释行，并手动指定列名 / Use skiprows and names to set column names
            df = pd.read_csv(
                self.assembly_summary_path,
                sep='\t',
                skiprows=header_line + 1,
                names=header_cols,
                low_memory=False
            )
            filtered_df = df[
                (df['organism_name'].str.contains(genus, case=False, na=False)) &
                (df['assembly_level'] == level)
            ]
            if len(filtered_df) == 0:
                print(f"⚠️  No genomes found for {genus} with level '{level}'")
                return None
            summary_columns = [
                'assembly_accession', 'bioproject', 'biosample', 'wgs_master', 'refseq_category',
                'taxid', 'species_taxid', 'organism_name', 'infraspecific_name', 'isolate',
                'version_status', 'assembly_level', 'release_type', 'genome_rep', 'seq_rel_date',
                'asm_name', 'submitter', 'gbrs_paired_asm', 'paired_asm_comp', 'ftp_path',
                'excluded_from_refseq', 'relation_to_type_material'
            ]
            available_columns = [col for col in summary_columns if col in filtered_df.columns]
            summary_df = filtered_df[available_columns].copy()
            summary_df.to_csv(summary_file, index=False)
            print(f"📊 Extracted summary for {len(summary_df)} {genus} genomes to {save_dir}")
            print(f"📁 Summary saved to: {summary_file}")
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
    
    def download_single_genome(self, accession, output_dir):
        """
        下载单个基因组到指定目录 / Download a single genome to the specified directory
        Args:
            accession (str): 基因组accession号 / Genome accession number
            output_dir (Path): 输出目录 / Output directory
        Returns:
            bool: 是否成功 / Success status
        """
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # 创建zip文件夹用于存放下载的zip文件 / Create zip folder to store downloaded zip files
        zip_folder = output_dir / "zip"
        zip_folder.mkdir(exist_ok=True)
        
        filename = str(zip_folder / f"{accession}.zip")
        try:
            print(f"🔄 Downloading genome {accession} to {filename} ...")
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
    
    def batch_download_genomes(self, listfile, output_dir, batch_size=10):
        """
        批量下载基因组到指定目录 / Download genomes in batches to the specified directory
        Args:
            listfile (str): accession号列表文件 / Path to file containing accession numbers
            output_dir (Path): 输出目录 / Output directory
            batch_size (int): 并行下载数 / Number of concurrent downloads
        Returns:
            bool: 是否成功 / Success status
        """
        with open(listfile, 'r') as f:
            accessions = [line.strip().split('\t')[0] for line in f.readlines()]
        total_accessions = len(accessions)
        print(f"📊 Total genomes to download: {total_accessions}")
        success_count = 0
        with ThreadPoolExecutor(max_workers=batch_size) as executor:
            futures = [executor.submit(self.download_single_genome, acc, output_dir) for acc in accessions]
            for future in as_completed(futures):
                if future.result():
                    success_count += 1
        print(f"📈 Download completed: {success_count}/{total_accessions} successful")
        return success_count > 0
    
    def unzip_files(self, genus_dir):
        """
        解压指定目录中的所有zip文件 / Extract all zip files in the specified directory
        
        Args:
            genus_dir (Path): 包含zip文件的目录 / Directory containing zip files
            
        Returns:
            bool: 成功状态 / Success status
        """
        print(f"🔄 [{genus_dir.name}] Extracting zip files...")
        
        # 从zip文件夹中获取zip文件 / Get zip files from zip folder
        zip_folder = genus_dir / "zip"
        if not zip_folder.exists():
            print(f"⚠️  [{genus_dir.name}] Zip folder not found")
            return False
        
        zip_files = list(zip_folder.glob("*.zip"))
        if not zip_files:
            print(f"⚠️  [{genus_dir.name}] No zip files found")
            return False
        
        for zip_file in zip_files:
            try:
                # 创建解压目录（移除.zip扩展名）/ Create extraction directory (remove .zip extension)
                extract_dir = genus_dir / zip_file.stem
                with zipfile.ZipFile(zip_file, 'r') as zip_ref:
                    zip_ref.extractall(extract_dir)
                print(f"✅ Extracted: {zip_file.name}")
            except Exception as e:
                print(f"❌ Failed to extract {zip_file.name}: {e}")
                return False
        
        return True
    
    def unzip_files_to_temp(self, download_dir, temp_extract_dir):
        """
        解压指定目录中的所有zip文件到临时目录 / Extract all zip files to temporary directory
        
        Args:
            download_dir (Path): 包含zip文件的下载目录 / Directory containing zip files
            temp_extract_dir (Path): 临时解压目录 / Temporary extraction directory
            
        Returns:
            bool: 成功状态 / Success status
        """
        print(f"🔄 [{download_dir.name}] Extracting zip files to temporary directory...")
        
        # 从zip文件夹中获取zip文件 / Get zip files from zip folder
        zip_folder = download_dir / "zip"
        if not zip_folder.exists():
            print(f"⚠️  [{download_dir.name}] Zip folder not found")
            return False
        
        zip_files = list(zip_folder.glob("*.zip"))
        if not zip_files:
            print(f"⚠️  [{download_dir.name}] No zip files found")
            return False
        
        # 确保临时解压目录存在 / Ensure temporary extraction directory exists
        temp_extract_dir.mkdir(parents=True, exist_ok=True)
        
        for zip_file in zip_files:
            try:
                # 解压到临时目录（移除.zip扩展名）/ Extract to temporary directory (remove .zip extension)
                extract_subdir = temp_extract_dir / zip_file.stem
                with zipfile.ZipFile(zip_file, 'r') as zip_ref:
                    zip_ref.extractall(extract_subdir)
                print(f"✅ Extracted to temp: {zip_file.name}")
            except Exception as e:
                print(f"❌ Failed to extract {zip_file.name}: {e}")
                return False
        
        return True
    
    def collect_fna_files(self, genus_dir):
        """
        收集所有fna文件到genus目录根目录 / Collect all fna files to the genus directory root
        
        Args:
            genus_dir (Path): 属目录 / Genus directory
            
        Returns:
            bool: 是否成功 / Success status
        """
        print(f"🔄 [{genus_dir.name}] Collecting fna files...")
        
        # 查找所有解压目录下的fna文件 / Find all fna files in extracted directories
        fna_files = list(genus_dir.rglob("*.fna"))
        
        if not fna_files:
            print(f"⚠️  [{genus_dir.name}] No fna files found")
            return False
        
        # 移动所有fna文件到genus目录根目录 / Move all fna files to genus directory root
        success_count = 0
        for fna_file in fna_files:
            if fna_file.parent != genus_dir:  # 只移动不在根目录的文件 / Only move files not in root
                dest_file = genus_dir / fna_file.name
                try:
                    shutil.move(str(fna_file), str(dest_file))
                    print(f"✅ Moved file: {fna_file.name}")
                    success_count += 1
                except Exception as e:
                    print(f"❌ Failed to move file {fna_file.name}: {e}")
        
        return success_count > 0
    
    def organize_genus_files(self, genus_dir, target_subdir, prefix):
        """
        组织genus文件到最终目录结构 / Organize genus files into final directory structure
        
        Args:
            genus_dir (Path): 临时genus目录 / Temporary genus directory
            target_subdir (Path): 目标子目录 / Target subdirectory
            prefix (str): 文件前缀 / File prefix
            
        Returns:
            bool: 是否成功 / Success status
        """
        print(f"🔄 Organizing files with prefix {prefix}_ to {target_subdir}")
        
        # 创建目标子目录 / Create target subdirectory
        target_subdir.mkdir(parents=True, exist_ok=True)
        
        # 获取所有fna文件 / Get all fna files
        fna_files = list(genus_dir.glob("*.fna"))
        if not fna_files:
            print(f"⚠️  No fna files found for organization")
            return False
        
        # 移动并重命名文件 / Move and rename files
        success_count = 0
        for fna_file in fna_files:
            original_name = fna_file.name
            new_name = f"{prefix}_{original_name}"
            target_file = target_subdir / new_name
            
            try:
                shutil.move(str(fna_file), str(target_file))
                print(f"✅ Organized file: {original_name} -> {new_name}")
                success_count += 1
            except Exception as e:
                print(f"❌ Failed to organize file {original_name}: {e}")
        
        print(f"✅ File organization completed: {success_count} files moved to {target_subdir}")
        return success_count > 0
    
    def process_genus_complete(self, genus, level="Complete Genome", is_outgroup=False, target_genus_safe=None):
        """
        完整处理一个genus：下载、解压、组织文件 / Complete processing of a genus: download, extract, organize files
        
        新的组织策略 / New organization strategy:
        - 目标菌: target_genus/data/target_genus_prefix/ / Target genus: target_genus/data/target_genus_prefix/
        - 外群: target_genus/data/outgroup/ / Outgroups: target_genus/data/outgroup/
        - 解压文件: target_genus/extracted/ / Extracted files: target_genus/extracted/
        - Summary文件: 每个genus都生成自己的summary / Summary files: each genus generates its own summary
        
        Args:
            genus (str): 属名 / Genus name
            level (str): 组装级别过滤 / Assembly level filter
            is_outgroup (bool): 是否为外群 / Whether this is an outgroup genus
            target_genus_safe (str): 目标属名（下划线安全格式），仅outgroup时用 / Target genus name (safe, with underscores), used for outgroup
        Returns:
            tuple: (genus_name, success_status)
        """
        genus_safe = genus.replace(" ", "_")
        print(f"\n🚀 Starting processing {'outgroup ' if is_outgroup else 'target '}genus: {genus}")
        
        # 确定工作目录和最终目标目录 / Determine working and final target directories
        if is_outgroup and target_genus_safe:
            # 外群：在目标菌目录内工作，最终移动到目标菌的outgroup子目录 / Outgroup: work in target dir, final move to target's outgroup subdir
            target_base_dir = self.work_dir / target_genus_safe
            temp_genus_dir = target_base_dir / "extracted" / genus_safe  # 解压到extracted子目录 / Extract to extracted subdir
            final_target_dir = target_base_dir / "data" / "outgroup"
            # 外群的summary文件也保存在目标菌目录中 / Outgroup summary files also saved in target genus directory
            summary_target_dir = target_base_dir
            prefix = "out"
        else:
            # 目标菌：在自己的目录下工作，最终存储到自己的data子目录 / Target genus: work in own dir, final storage in own data subdir
            target_base_dir = self.work_dir / genus_safe
            temp_genus_dir = target_base_dir / "extracted" / genus_safe  # 解压到extracted子目录 / Extract to extracted subdir
            final_target_dir = target_base_dir / "data" / genus_safe[:3].capitalize()
            summary_target_dir = target_base_dir
            prefix = genus_safe[:3].capitalize()
        
        print(f"📂 Target base directory: {target_base_dir}")
        print(f"📂 Temporary extraction directory: {temp_genus_dir}")
        print(f"📁 Final target directory: {final_target_dir}")
        print(f"📊 Summary target directory: {summary_target_dir}")
        print(f"🏷️  File prefix: {prefix}_")
        
        # 创建必要的目录 / Create necessary directories
        target_base_dir.mkdir(parents=True, exist_ok=True)
        temp_genus_dir.mkdir(parents=True, exist_ok=True)
        
        # 下载基因组到目标基础目录 / Download genomes to target base directory
        if not self.download_genus_with_summary(genus, level, target_base_dir, summary_target_dir):
            return genus, False
        
        # 解压文件到临时目录 / Extract files to temporary directory
        if not self.unzip_files_to_temp(target_base_dir, temp_genus_dir):
            return genus, False
        
        # 收集fna文件 / Collect fna files
        if not self.collect_fna_files(temp_genus_dir):
            return genus, False
        
        # 组织文件到最终目录 / Organize files to final directory
        if not self.organize_genus_files(temp_genus_dir, final_target_dir, prefix):
            return genus, False
        
        # 清理临时解压目录 / Clean up temporary extraction directory
        try:
            if temp_genus_dir.exists():
                shutil.rmtree(temp_genus_dir)
                print(f"🧹 Cleaned up temporary extraction directory: {temp_genus_dir}")
        except Exception as e:
            print(f"⚠️  Warning: Failed to clean up temporary directory {temp_genus_dir}: {e}")
        
        print(f"🎉 [{genus}] Complete processing finished!")
        return genus, True
    
    def download_genus_with_summary(self, genus, level="Complete Genome", output_dir=None, summary_dir=None):
        """
        下载指定genus的基因组并提取summary信息到指定目录 / Download genomes for a specific genus and extract summary to specified directory
        Args:
            genus (str): 目标属名 / Target genus name
            level (str): 组装级别过滤 / Assembly level filter
            output_dir (Path): 输出目录 / Output directory
            summary_dir (Path): summary文件保存目录 / Summary file save directory
        Returns:
            bool: 是否成功 / Success status
        """
        genus_safe = genus.replace(" ", "_")
        if output_dir is not None:
            genus_dir = output_dir
        else:
            genus_dir = self.work_dir / genus_safe
        genus_dir.mkdir(exist_ok=True)
        
        # 确定summary文件保存目录 / Determine summary file save directory
        if summary_dir is not None:
            summary_save_dir = summary_dir
        else:
            summary_save_dir = genus_dir
        summary_save_dir.mkdir(exist_ok=True)
        
        log_file = genus_dir / "download.log"
        not_found_log = self.work_dir / "not_found.log"
        
        # 提取基因组列表 / Extract genome list
        list_file = self.extract_genomes_by_level(genus, level)
        if not list_file:
            with open(not_found_log, "a") as nf:
                nf.write(f"[{datetime.now()}] NOT FOUND: {genus} (level: {level})\n")
            print(f"❌ No genomes found for {genus} with level '{level}'")
            return False
        
        # 提取基因组summary信息到指定目录 / Extract genome summary information to specified directory
        print(f"📊 Extracting genome summary information for {genus} to {summary_save_dir}...")
        summary_file = self.extract_genome_summary_to_dir(genus, level, summary_save_dir)
        if summary_file:
            print(f"✅ Genome summary extracted successfully to {summary_file}")
        else:
            print(f"⚠️  Failed to extract genome summary, but continuing with download")
        
        try:
            # 下载基因组 / Download genomes
            success = self.batch_download_genomes(list_file, genus_dir)
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
    
    def download_genus(self, genus, level="Complete Genome", output_dir=None):
        """
        下载指定genus的基因组并提取summary信息 / Download genomes for a specific genus and extract summary information
        Args:
            genus (str): 目标属名 / Target genus name
            level (str): 组装级别过滤 / Assembly level filter
            output_dir (Path): 输出目录 / Output directory
        Returns:
            bool: 是否成功 / Success status
        """
        genus_safe = genus.replace(" ", "_")
        if output_dir is not None:
            genus_dir = output_dir
        else:
            genus_dir = self.work_dir / genus_safe
        genus_dir.mkdir(exist_ok=True)
        
        log_file = genus_dir / "download.log"
        not_found_log = self.work_dir / "not_found.log"
        
        # 提取基因组列表 / Extract genome list
        list_file = self.extract_genomes_by_level(genus, level)
        if not list_file:
            with open(not_found_log, "a") as nf:
                nf.write(f"[{datetime.now()}] NOT FOUND: {genus} (level: {level})\n")
            print(f"❌ No genomes found for {genus} with level '{level}'")
            return False
        
        # 提取基因组summary信息 / Extract genome summary information
        print(f"📊 Extracting genome summary information for {genus}...")
        summary_file = self.extract_genome_summary(genus, level)
        if summary_file:
            print(f"✅ Genome summary extracted successfully")
        else:
            print(f"⚠️  Failed to extract genome summary, but continuing with download")
        
        try:
            # 下载基因组 / Download genomes
            success = self.batch_download_genomes(list_file, genus_dir)
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
        下载目标genus和外群 / Download target genus and outgroup genomes
        
        新的文件组织策略 / New file organization strategy:
        target_genus/
        ├── data/
        │   ├── Target_genus_prefix/    # 目标菌基因组文件 / Target genus genome files
        │   │   ├── Tar_genome1.fna
        │   │   └── Tar_genome2.fna
        │   └── outgroup/               # 外群基因组文件 / Outgroup genome files
        │       ├── out_genome1.fna
        │       └── out_genome2.fna
        ├── zip/                        # 下载的zip文件 / Downloaded zip files
        │   ├── target_genome1.zip
        │   └── outgroup_genome1.zip
        ├── target_genus_genome_summary.csv      # 目标菌群summary / Target genus summary
        ├── outgroup1_genus_genome_summary.csv   # 外群1 summary / Outgroup1 summary
        └── outgroup2_genus_genome_summary.csv   # 外群2 summary / Outgroup2 summary
        
        注意：extracted/ 目录用于临时解压，处理完成后会被清理
        Note: extracted/ directory is used for temporary extraction and cleaned up after processing
        
        Args:
            target_genus (str): 目标属名 / Target genus name
            outgroup_genera (list): 外群属名列表 / List of outgroup genus names
            level (str): 组装级别过滤 / Assembly level filter
            threads (int): 并行线程数 / Number of parallel threads
        Returns:
            bool: 是否成功 / Success status
        """
        print(f"⚙️ Starting genome download with improved file organization")
        print(f"📋 Target genus: {target_genus}")
        if outgroup_genera:
            print(f"🔗 Outgroups: {', '.join(outgroup_genera)}")
        
        target_genus_safe = target_genus.replace(" ", "_")
        
        print(f"\n📁 Final directory structure will be:")
        print(f"   {target_genus_safe}/")
        print(f"   ├── data/")
        print(f"   │   ├── {target_genus_safe[:3].capitalize()}/    # Target genus genomes")
        print(f"   │   └── outgroup/                                # Outgroup genomes")
        print(f"   ├── zip/                                         # Downloaded zip files")
        print(f"   ├── {target_genus_safe}_genome_summary.csv       # Target genus summary")
        if outgroup_genera:
            for outgroup in outgroup_genera:
                outgroup_safe = outgroup.replace(" ", "_")
                print(f"   ├── {outgroup_safe}_genome_summary.csv       # Outgroup summary")
        print(f"   ")
        print(f"   Note: extracted/ directory is used for temporary extraction and cleaned up")
        
        # 首先处理目标菌群 / First process target genus
        print(f"\n🔄 Step 1: Processing target genus {target_genus}...")
        target_result = self.process_genus_complete(target_genus, level, False, None)
        target_success = target_result[1]
        
        if not target_success:
            print(f"❌ Target genus {target_genus} processing failed")
            return False
        
        # 然后处理外群 / Then process outgroups
        outgroup_success_count = 0
        total_outgroups = len(outgroup_genera) if outgroup_genera else 0
        
        if outgroup_genera:
            print(f"\n🔄 Step 2: Processing {total_outgroups} outgroups...")
            for i, outgroup in enumerate(outgroup_genera, 1):
                print(f"\n   🔄 Processing outgroup {i}/{total_outgroups}: {outgroup}...")
                genus_result, success = self.process_genus_complete(outgroup, level, True, target_genus_safe)
                if success:
                    outgroup_success_count += 1
                    print(f"   ✅ Outgroup {outgroup} processed successfully")
                else:
                    print(f"   ❌ Outgroup {outgroup} processing failed")
        
        # 验证最终文件组织 / Verify final file organization
        self._verify_file_organization(target_genus_safe, outgroup_genera)
        
        print(f"\n🏁 Download and organization completed!")
        print(f"📊 Results summary:")
        print(f"   ✅ Target genus: {1 if target_success else 0}/1")
        print(f"   ✅ Outgroups: {outgroup_success_count}/{total_outgroups}")
        
        success_count = (1 if target_success else 0) + outgroup_success_count
        return success_count > 0
    
    def _verify_file_organization(self, target_genus_safe, outgroup_genera=None):
        """
        验证最终文件组织结构 / Verify final file organization structure
        
        Args:
            target_genus_safe (str): 目标属名（下划线安全格式）/ Target genus name (safe, with underscores)
            outgroup_genera (list): 外群属名列表 / List of outgroup genus names
        """
        print(f"\n📋 Verifying file organization for {target_genus_safe}...")
        
        target_base_dir = self.work_dir / target_genus_safe
        target_data_dir = target_base_dir / "data"
        
        if not target_data_dir.exists():
            print(f"   ⚠️  Data directory not found: {target_data_dir}")
            return
        
        # 检查目标菌文件 / Check target genus files
        target_genus_dir = target_data_dir / target_genus_safe[:3].capitalize()
        if target_genus_dir.exists():
            target_files = list(target_genus_dir.glob("*.fna"))
            print(f"   📂 Target genus files ({len(target_files)} files):")
            for file in target_files[:3]:  # 显示前3个文件 / Show first 3 files
                print(f"      ✓ {file.name}")
            if len(target_files) > 3:
                print(f"      ... and {len(target_files) - 3} more files")
        else:
            print(f"   ⚠️  Target genus directory not found: {target_genus_dir}")
        
        # 检查外群文件 / Check outgroup files
        outgroup_dir = target_data_dir / "outgroup"
        if outgroup_dir.exists():
            outgroup_files = list(outgroup_dir.glob("*.fna"))
            print(f"   📂 Outgroup files ({len(outgroup_files)} files):")
            for file in outgroup_files[:3]:  # 显示前3个文件 / Show first 3 files
                print(f"      ✓ {file.name}")
            if len(outgroup_files) > 3:
                print(f"      ... and {len(outgroup_files) - 3} more files")
        else:
            print(f"   📂 Outgroup directory not found (no outgroups processed)")
        
        # 检查zip目录 / Check zip directory
        zip_dir = target_base_dir / "zip"
        if zip_dir.exists():
            zip_files = list(zip_dir.glob("*.zip"))
            print(f"   ✅ zip directory exists: {len(zip_files)} files")
        else:
            print(f"   ⚠️  zip directory not found")
        
        # 检查目标菌summary文件 / Check target genus summary file
        target_summary_file = target_base_dir / f"{target_genus_safe}_genome_summary.csv"
        if target_summary_file.exists():
            print(f"   ✅ Target genus summary file: {target_summary_file.name}")
        else:
            print(f"   ⚠️  Target genus summary file not found: {target_summary_file.name}")
        
        # 检查外群summary文件 / Check outgroup summary files
        if outgroup_genera:
            print(f"   📊 Checking outgroup summary files:")
            for outgroup in outgroup_genera:
                outgroup_safe = outgroup.replace(" ", "_")
                outgroup_summary_file = target_base_dir / f"{outgroup_safe}_genome_summary.csv"
                if outgroup_summary_file.exists():
                    print(f"      ✅ {outgroup_safe}_genome_summary.csv")
                else:
                    print(f"      ⚠️  {outgroup_safe}_genome_summary.csv not found")
        else:
            print(f"   📊 No outgroups to check summary files for")


def main():
    """
    基因组下载器的命令行接口 / Command line interface for GenomeDownloader
    """
    import argparse
    
    parser = argparse.ArgumentParser(description="Download bacterial genomes using assembly_summary.txt")
    parser.add_argument("genus", help="Target genus name")
    parser.add_argument("--outgroup", nargs="*", help="Outgroup genus names")
    parser.add_argument("--level", default="Complete Genome", 
                       help="Assembly level (default: 'Complete Genome')")
    parser.add_argument("--threads", type=int, default=64, 
                       help="Number of threads (default: 4)")
    parser.add_argument("--work-dir", help="Working directory")
    parser.add_argument("--assembly-summary", help="Path to assembly_summary.txt file")
    
    args = parser.parse_args()
    
    try:
        # 创建下载器 / Create downloader
        downloader = GenomeDownloader(args.work_dir, args.assembly_summary)
        
        # 下载基因组 / Download genomes
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
