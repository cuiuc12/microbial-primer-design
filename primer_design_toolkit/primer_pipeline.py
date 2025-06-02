#!/usr/bin/env python3
"""
Complete Bacterial Genome Analysis Pipeline: Download -> Prokka -> Roary -> Primer Design
Supports resuming from any step in the workflow
"""

import os
import sys
import subprocess
import glob
import shutil
import argparse
from datetime import datetime
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

class GenomePipeline:
    """
    Complete pipeline for bacterial genome analysis and primer design
    
    Workflow:
    1. Download genomes from NCBI
    2. Annotate with Prokka 
    3. Pangenome analysis with Roary
    4. Identify species-specific genes
    5. Extract gene sequences
    6. Multiple sequence alignment with MAFFT
    7. Primer design with Primer3
    8. Parse Primer3 results
    9. Quality scoring and ranking
    """
    def __init__(self, genus, outgroup_genera=None, level="Complete Genome", threads=4, fast_mode=True, assembly_summary_path=None, max_genes=None):
        """
        Initialize the pipeline
        
        Args:
            genus (str): Target genus name
            outgroup_genera (list): List of outgroup genus names
            level (str): Assembly level filter
            threads (int): Number of parallel threads
            fast_mode (bool): Use fast mode for Roary analysis
            assembly_summary_path (str): Path to assembly_summary.txt file
            max_genes (int): Maximum number of genes to process (None for all genes)
        """
        self.genus = genus
        self.genus_safe = genus.replace(" ", "_")
        self.outgroup_genera = outgroup_genera or []
        self.outgroup_genera_safe = [g.replace(" ", "_") for g in (outgroup_genera or [])]
        self.level = level
        self.threads = threads
        self.fast_mode = fast_mode  # Control roary speed mode
        self.assembly_summary_path = assembly_summary_path
        self.max_genes = max_genes  # Maximum number of genes to process
        self.work_dir = Path.cwd()
        self.genus_dir = self.work_dir / self.genus_safe
        self.data_dir = self.genus_dir / "data"
        self.target_dir = self.data_dir / self.genus_safe[:3].capitalize()
        self.outgroup_dir = self.data_dir / "outgroup"
        self.prokka_dir = self.genus_dir / "prokka_results"
        self.roary_dir = self.genus_dir / "roary_results"
        
    def log(self, message, level="INFO"):
        """
        Log messages with timestamp
        
        Args:
            message (str): Log message
            level (str): Log level (INFO, ERROR, etc.)
        """
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        print(f"[{timestamp}] {level}: {message}")
        
        # 确保日志目录存在
        self.genus_dir.mkdir(parents=True, exist_ok=True)
        
        # Write to log file
        log_file = self.genus_dir / "pipeline.log"
        with open(log_file, "a", encoding="utf-8") as f:
            f.write(f"[{timestamp}] {level}: {message}\n")
    
    def check_conda_env(self, env_name):
        """检查conda环境是否存在"""
        try:
            result = subprocess.run(
                ["conda", "env", "list"], 
                capture_output=True, 
                text=True, 
                check=True
            )
            return env_name in result.stdout
        except subprocess.CalledProcessError:
            return False
    
    def activate_conda_env(self, env_name):
        """激活conda环境 / Activate conda environment"""
        if not self.check_conda_env(env_name):
            self.log(f"❌ Conda environment '{env_name}' not found", "ERROR")
            return False
        
        self.log(f"🔄 Activating conda environment: {env_name}")
        return True
    
    def run_command_in_env(self, command, env_name, cwd=None):
        """在指定conda环境中运行命令 / Run command in specified conda environment"""
        if not self.activate_conda_env(env_name):
            return False
            
        # 构建在conda环境中运行的命令 / Build command to run in conda environment
        if isinstance(command, list):
            cmd_str = " ".join(command)
        else:
            cmd_str = command
            
        full_command = f"conda run -n {env_name} {cmd_str}"
        
        try:
            self.log(f"🔄 Executing command: {cmd_str}")
            result = subprocess.run(
                full_command,
                shell=True,
                cwd=cwd or self.work_dir,
                check=True,
                capture_output=False
            )
            return True
        except subprocess.CalledProcessError as e:
            self.log(f"❌ Command execution failed: {e}", "ERROR")
            return False
    
    def step1_download(self):
        """步骤1: 下载基因组 / Step 1: Download genomes"""
        self.log("🚀 Starting Step 1: Genome download")
        
        # 检查是否已经下载完成 / Check if download is already completed
        if self.target_dir.exists() and list(self.target_dir.glob("*.fna")):
            self.log("✅ Target genomes already exist, skipping download")
            target_downloaded = True
        else:
            target_downloaded = False
            
        if self.outgroup_genera and self.outgroup_dir.exists() and list(self.outgroup_dir.glob("*.fna")):
            self.log("✅ Outgroup genomes already exist, skipping download")
            outgroup_downloaded = True
        else:
            outgroup_downloaded = False
            
        if target_downloaded and (not self.outgroup_genera or outgroup_downloaded):
            self.log("✅ All genomes already downloaded")
            # 检查并显示已有的summary文件 / Check and display existing summary files
            self._check_existing_summary_files()
            return True
        
        # 使用新的GenomeDownloader / Use new GenomeDownloader
        try:
            from .genome_downloader import GenomeDownloader
            
            downloader = GenomeDownloader(
                work_dir=str(self.work_dir),
                assembly_summary_path=self.assembly_summary_path
            )
            
            # 下载目标基因组和外群 / Download target genomes and outgroups
            success = downloader.download_with_outgroup(
                target_genus=self.genus,
                outgroup_genera=self.outgroup_genera,
                level=self.level,
                threads=self.threads
            )
            
            if success:
                self.log("✅ Step 1 completed: Genome download")
                # 显示summary文件位置 / Display summary file locations
                self._display_summary_files()
                return True
            else:
                self.log("❌ Genome download failed", "ERROR")
                return False
                
        except Exception as e:
            self.log(f"❌ Error during download: {e}", "ERROR")
            return False
    
    def _check_existing_summary_files(self):
        """检查并显示已存在的summary文件 / Check and display existing summary files"""
        summary_files = []
        
        # Check target genus summary file
        target_summary = self.genus_dir / f"{self.genus_safe}_genome_summary.csv"
        if target_summary.exists():
            summary_files.append(f"Target genome summary: {target_summary}")
        
        # Check outgroup summary files
        for outgroup, outgroup_safe in zip(self.outgroup_genera, self.outgroup_genera_safe):
            outgroup_dir = self.work_dir / outgroup_safe
            outgroup_summary = outgroup_dir / f"{outgroup_safe}_genome_summary.csv"
            if outgroup_summary.exists():
                summary_files.append(f"Outgroup genome summary: {outgroup_summary}")
        
        if summary_files:
            self.log("📊 Genome summary files saved:")
            for summary_file in summary_files:
                self.log(f"   {summary_file}")
        else:
            self.log("⚠️  No genome summary file found")
    
    def _display_summary_files(self):
        """显示新生成的summary文件位置 / Display newly generated summary file locations"""
        summary_files = []
        
        # 检查目标genus的summary文件 / Check target genus summary file
        target_summary = self.genus_dir / f"{self.genus_safe}_genome_summary.csv"
        if target_summary.exists():
            summary_files.append(f"Target genome summary: {target_summary}")
        
        # 检查外群的summary文件 / Check outgroup summary files
        for outgroup, outgroup_safe in zip(self.outgroup_genera, self.outgroup_genera_safe):
            outgroup_dir = self.work_dir / outgroup_safe
            outgroup_summary = outgroup_dir / f"{outgroup_safe}_genome_summary.csv"
            if outgroup_summary.exists():
                summary_files.append(f"Outgroup genome summary: {outgroup_summary}")
        
        if summary_files:
            self.log("📊 Genome summary information saved:")
            for summary_file in summary_files:
                self.log(f"   {summary_file}")
        else:
            self.log("⚠️  No genome summary file found")
    
    def step2_prokka(self):
        """步骤2: Prokka注释 (并行优化版本) / Step 2: Prokka annotation (parallel optimized version)"""
        self.log("🚀 Starting Step 2: Prokka annotation")
        
        # 创建prokka结果目录 / Create prokka results directory
        self.prokka_dir.mkdir(exist_ok=True)
        
        # Collect all fna files that need annotation
        fna_files = []
        
        # 目标基因组文件 / Target genome files
        if self.target_dir.exists():
            fna_files.extend(list(self.target_dir.glob("*.fna")))
        
        # 外群文件 / Outgroup files
        if self.outgroup_dir.exists():
            fna_files.extend(list(self.outgroup_dir.glob("*.fna")))
        
        if not fna_files:
            self.log("❌ No fna files found for annotation", "ERROR")
            return False
        
        self.log(f"📋 Found {len(fna_files)} genomes for annotation")
        
        # 检查已完成的注释 / Check completed annotations
        remaining_files = []
        for fna_file in fna_files:
            sample_name = fna_file.stem
            gff_file = self.prokka_dir / sample_name / f"{sample_name}.gff"
            if gff_file.exists():
                self.log(f"✅ {sample_name} annotation already exists, skipping")
            else:
                remaining_files.append(fna_file)
        
        if not remaining_files:
            self.log("✅ All genome annotations completed")
            return True
        
        self.log(f"🔄 Need to annotate {len(remaining_files)} genomes")
        
        # 一次性检查prokka环境 / Check prokka environment once
        if not self.check_conda_env("prokka"):
            self.log("❌ Conda environment 'prokka' not found", "ERROR")
            return False
        
        # 优化并行策略和CPU分配 / Optimized parallel strategy and CPU allocation
        if len(remaining_files) <= 3:
            # 少量基因组：串行处理，每个使用更多CPU / Few genomes: serial processing with more CPUs each
            max_workers = 1
            cpus_per_prokka = min(self.threads, 16)  # 单个Prokka进程最多使用16个CPU / Single Prokka process uses max 16 CPUs
            self.log(f"🚀 Using serial processing strategy: {cpus_per_prokka} CPUs per genome")
        else:
            # 多个基因组：并行处理，但限制并发数以避免资源竞争
            # Multiple genomes: parallel processing with limited concurrency to avoid resource contention
            max_workers = min(max(1, self.threads // 8), len(remaining_files))  # 每个进程至少8个CPU / At least 8 CPUs per process
            cpus_per_prokka = max(1, self.threads // max_workers)
            cpus_per_prokka = min(cpus_per_prokka, 16)  # 限制最大CPU数 / Limit max CPU count
            self.log(f"🚀 Using parallel processing strategy: {max_workers} parallel processes, {cpus_per_prokka} CPUs each")
        
        # 并行执行prokka注释 / Execute prokka annotation in parallel
        def run_single_prokka(fna_file):
            """运行单个基因组的prokka注释 / Run prokka annotation for single genome"""
            sample_name = fna_file.stem
            output_dir = self.prokka_dir / sample_name
            
            cmd = [
                "conda", "run", "-n", "prokka",
                "prokka",
                "--outdir", str(output_dir),
                "--prefix", sample_name,
                "--kingdom", "Bacteria",
                "--cpus", str(cpus_per_prokka),  # 使用预计算的CPU数量 / Use pre-calculated CPU count
                "--force",
                str(fna_file)
            ]
            
            try:
                self.log(f"🔄 Starting annotation: {sample_name} (using {cpus_per_prokka} CPUs)")
                result = subprocess.run(
                    cmd,
                    cwd=self.work_dir,
                    check=True,
                    capture_output=True,
                    text=True
                )
                self.log(f"✅ Annotation completed: {sample_name}")
                return sample_name, True, None
            except subprocess.CalledProcessError as e:
                error_msg = f"stderr: {e.stderr[:200]}..." if e.stderr else str(e)
                self.log(f"❌ Annotation failed {sample_name}: {error_msg}", "ERROR")
                return sample_name, False, error_msg
        
        # 计算合适的并行数量 / Calculate appropriate parallel count
        # 总CPU数除以每个prokka进程的CPU数，但不超过基因组数量 / Total CPUs divided by CPUs per prokka process, but not exceeding genome count
        max_workers = min(self.threads, len(remaining_files))
        self.log(f"🚀 Using {max_workers} parallel processes for annotation")
        
        # 执行并行注释 / Execute parallel annotation
        success_count = 0
        failed_samples = []
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # 提交所有任务 / Submit all tasks
            future_to_file = {
                executor.submit(run_single_prokka, fna_file): fna_file 
                for fna_file in remaining_files
            }
            
            # 处理完成的任务 / Process completed tasks
            for future in as_completed(future_to_file):
                fna_file = future_to_file[future]
                try:
                    sample_name, success, error = future.result()
                    if success:
                        success_count += 1
                    else:
                        failed_samples.append(sample_name)
                except Exception as exc:
                    sample_name = fna_file.stem
                    self.log(f"❌ {sample_name} execution error: {exc}", "ERROR")
                    failed_samples.append(sample_name)
        
        # 报告结果 / Report results
        self.log(f"📊 Annotation completion statistics: {success_count}/{len(remaining_files)} successful")
        if failed_samples:
            self.log(f"❌ Failed samples: {', '.join(failed_samples)}", "ERROR")
            return False
        
        self.log("✅ Step 2 completed: Prokka annotation")
        return True
    
    def step3_roary(self):
        """步骤3: Roary泛基因组分析 (优化版本)"""
        self.log("🚀 开始步骤3: Roary泛基因组分析")
        
        # 检查是否已经完成 - 检查多个可能的结果目录
        result_files_to_check = [
            self.roary_dir / "gene_presence_absence.csv",
            self.work_dir / "roary_out" / "gene_presence_absence.csv"
        ]
        
        # 检查是否有现有的roary结果目录（包括带时间戳的）
        existing_roary_dirs = list(self.work_dir.glob("roary_*"))
        for roary_dir in existing_roary_dirs:
            if (roary_dir / "gene_presence_absence.csv").exists():
                self.log(f"✅ 发现已完成的Roary分析: {roary_dir}")
                return True
        
        for result_file in result_files_to_check:
            if result_file.exists():
                self.log(f"✅ Roary分析已完成: {result_file.parent}")
                return True
        
        # 收集所有gff文件
        gff_files = list(self.prokka_dir.glob("*/*.gff"))
        
        if not gff_files:
            self.log("❌ 没有找到gff文件，请先运行Prokka注释", "ERROR")
            return False
        
        self.log(f"📋 找到 {len(gff_files)} 个gff文件用于泛基因组分析")
        
        # 一次性检查roary环境
        if not self.check_conda_env("roary"):
            self.log("❌ Conda环境 'roary' 不存在", "ERROR")
            return False
        
        # 清理可能存在的空目录，避免roary创建时间戳目录
        if self.roary_dir.exists():
            try:
                # 如果目录为空或不包含主要结果文件，删除它
                if not any(self.roary_dir.iterdir()) or not (self.roary_dir / "gene_presence_absence.csv").exists():
                    shutil.rmtree(self.roary_dir)
                    self.log(f"🧹 清理空的roary目录: {self.roary_dir}")
            except Exception as e:
                self.log(f"⚠️  清理目录时出现警告: {e}")
        
        # 不要预先创建输出目录，让Roary自己创建
        # 移除这行: self.roary_dir.mkdir(exist_ok=True)
        
        # 运行Roary - 根据模式选择参数
        if self.fast_mode:
            # 快速模式：适合初步分析和大规模数据
            mode_desc = "快速模式 (-n)"
            mode_params = [
                "-e",           # 创建多序列比对
                "-n",           # 快速核心基因比对
            ]
        else:
            # 高质量模式：适合最终分析和引物设计
            mode_desc = "高质量模式 (--mafft)"
            mode_params = [
                "-e",           # 创建多序列比对
                "--mafft",      # 使用MAFFT进行高质量比对
            ]
        
        cmd = [
            "conda", "run", "-n", "roary",
            "roary",
        ] + mode_params + [
            "-p", str(self.threads),  # 线程数
            "-i", "95",     # 身份阈值95% (引物设计需要更高保守性)
            "-cd", "100",   # 核心基因定义阈值100% (引物目标必须在所有基因组中存在)
            "-v",           # 详细输出
            "-f", str(self.roary_dir)  # 输出目录 - Roary会创建这个目录
        ] + [str(f) for f in gff_files]
        
        try:
            self.log(f"🔄 运行Roary泛基因组分析 - {mode_desc}")
            self.log(f"📁 Roary将创建输出目录: {self.roary_dir}")
            result = subprocess.run(
                cmd,
                cwd=self.work_dir,
                check=True,
                capture_output=False  # 显示实时输出
            )
            
            # 验证结果文件是否成功创建
            if (self.roary_dir / "gene_presence_absence.csv").exists():
                self.log("✅ 步骤3完成: Roary泛基因组分析")
                return True
            else:
                self.log("❌ Roary运行完成但未找到预期的结果文件", "ERROR")
                return False
                
        except subprocess.CalledProcessError as e:
            self.log(f"❌ Roary分析失败: {e}", "ERROR")
            return False
    
    def step4_find_specific_genes(self):
        """步骤4: 找出只在目标genus中存在而在外群中不存在的基因 (并行优化版本)"""
        self.log("🚀 开始步骤4: 筛选目标特异性基因")
        
        # 检查roary结果是否存在
        gene_presence_file = None
        roary_dirs = [self.roary_dir]
        
        # 检查是否有其他roary结果目录
        existing_roary_dirs = list(self.work_dir.glob("roary_*"))
        roary_dirs.extend(existing_roary_dirs)
        
        for roary_dir in roary_dirs:
            potential_file = roary_dir / "gene_presence_absence.csv"
            if potential_file.exists():
                gene_presence_file = potential_file
                break
        
        if not gene_presence_file:
            self.log("❌ 没有找到gene_presence_absence.csv文件，请先运行Roary分析", "ERROR")
            return False
        
        self.log(f"📋 使用基因存在/缺失矩阵: {gene_presence_file}")
        
        # 分析基因存在/缺失模式
        import pandas as pd
        df = pd.read_csv(gene_presence_file, low_memory=False)
        
        # 获取所有列名（样本名）
        sample_columns = [col for col in df.columns if col not in 
                         ['Gene', 'Non-unique Gene name', 'Annotation', 'No. isolates', 'No. sequences', 
                          'Avg sequences per isolate', 'Genome Fragment', 'Order within Fragment', 
                          'Accessory Fragment', 'Accessory Order with Fragment', 'QC', 'Min group size nuc', 
                          'Max group size nuc', 'Avg group size nuc']]
        
        self.log(f"📊 发现 {len(sample_columns)} 个样本")
        
        # 简单而可靠的样本分类：目标 vs 外群
        target_samples = []
        outgroup_samples = []
        
        # 目标genus的前缀
        target_prefix = self.genus[:3].capitalize() + "_"
        
        self.log(f"🎯 目标前缀: {target_prefix}")
        
        # 分类样本：是目标前缀的为目标，其他都是外群
        for sample in sample_columns:
            if sample.startswith(target_prefix):
                target_samples.append(sample)
            else:
                outgroup_samples.append(sample)
        
        # 日志输出识别结果
        self.log(f"🎯 目标基因组: {len(target_samples)} 个")
        if target_samples:
            self.log(f"   样本: {', '.join(target_samples[:3])}{'...' if len(target_samples) > 3 else ''}")
        
        self.log(f"🔗 外群基因组: {len(outgroup_samples)} 个")
        if outgroup_samples:
            self.log(f"   样本: {', '.join(outgroup_samples[:3])}{'...' if len(outgroup_samples) > 3 else ''}")
        
        if not target_samples:
            self.log("❌ 没有找到目标基因组样本", "ERROR")
            return False
        
        # 最终统计
        self.log(f"📊 最终分类结果:")
        self.log(f"   🎯 目标基因组: {len(target_samples)} 个")
        self.log(f"   🔗 外群基因组: {len(outgroup_samples)} 个")
        
        # 并行筛选特异性基因
        def check_gene_specificity(gene_data):
            """并行检查单个基因的特异性"""
            try:
                idx, row = gene_data
                gene_name = row['Gene']
                annotation = str(row.get('Annotation', '')).lower()
                
                # 检查是否在所有目标基因组中存在
                target_present = all(pd.notna(row[sample]) and str(row[sample]).strip() != '' 
                                   for sample in target_samples)
                
                # 检查是否在所有外群中不存在
                if outgroup_samples:
                    outgroup_absent = all(pd.isna(row[sample]) or str(row[sample]).strip() == '' 
                                        for sample in outgroup_samples)
                else:
                    outgroup_absent = True  # 如果没有外群，则认为符合条件
                
                is_specific = target_present and outgroup_absent
                is_hypothetical = 'hypothetical protein' in annotation
                
                return idx, {
                    'gene_name': gene_name,
                    'annotation': annotation,
                    'is_specific': is_specific,
                    'is_hypothetical': is_hypothetical,
                    'target_present': target_present,
                    'outgroup_absent': outgroup_absent
                }, True, None
                
            except Exception as e:
                return idx, {}, False, str(e)
        
        # 并行处理所有基因
        self.log("🔄 开始并行筛选基因...")
        max_workers = min(self.threads, len(df))
        self.log(f"🚀 使用 {max_workers} 个并行进程进行基因筛选")
        
        gene_results = {}
        failed_count = 0
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # 提交所有基因检查任务
            future_to_gene = {
                executor.submit(check_gene_specificity, (idx, row)): idx 
                for idx, row in df.iterrows()
            }
            
            # 处理完成的任务
            for future in as_completed(future_to_gene):
                gene_idx = future_to_gene[future]
                try:
                    idx, gene_info, success, error = future.result()
                    if success:
                        gene_results[idx] = gene_info
                    else:
                        failed_count += 1
                        if error:
                            self.log(f"❌ 基因 {idx} 筛选失败: {error}", "ERROR")
                except Exception as exc:
                    failed_count += 1
                    self.log(f"❌ 基因 {gene_idx} 筛选出错: {exc}", "ERROR")
        
        # 收集筛选结果
        all_specific_genes = []
        hypothetical_specific_genes = []
        
        for idx, gene_info in gene_results.items():
            if gene_info.get('is_specific', False):
                gene_name = gene_info['gene_name']
                all_specific_genes.append(gene_name)
                
                if gene_info.get('is_hypothetical', False):
                    hypothetical_specific_genes.append(gene_name)
        
        self.log(f"✅ 发现 {len(all_specific_genes)} 个目标特异性基因")
        self.log(f"🎯 其中 {len(hypothetical_specific_genes)} 个为hypothetical protein")
        if failed_count > 0:
            self.log(f"❌ 失败基因数: {failed_count}")
        
        # 优先使用hypothetical protein，如果数量太少则使用所有特异性基因
        min_genes_threshold = 5  # 至少需要5个基因才认为足够
        
        if len(hypothetical_specific_genes) >= min_genes_threshold:
            selected_genes = hypothetical_specific_genes
            self.log(f"🔍 优先使用 {len(selected_genes)} 个hypothetical protein基因进行引物设计")
        else:
            selected_genes = all_specific_genes
            self.log(f"⚠️  hypothetical protein基因数量不足({len(hypothetical_specific_genes)} < {min_genes_threshold})")
            self.log(f"🔄 使用所有 {len(selected_genes)} 个特异性基因进行引物设计")
        
        # 保存选定的特异性基因列表
        specific_genes_file = self.genus_dir / "specific_genes.txt"
        with open(specific_genes_file, 'w') as f:
            for gene in selected_genes:
                f.write(f"{gene}\n")
        
        # 同时保存详细的基因信息
        detailed_genes_file = self.genus_dir / "specific_genes_detailed.csv"
        selected_df = df[df['Gene'].isin(selected_genes)][['Gene', 'Annotation']].copy()
        selected_df.to_csv(detailed_genes_file, index=False)
        
        self.log(f"📁 选定基因列表保存到: {specific_genes_file}")
        self.log(f"📁 详细基因信息保存到: {detailed_genes_file}")
        
        if len(selected_genes) == 0:
            self.log("⚠️  没有找到目标特异性基因，无法进行引物设计", "ERROR")
            return False
        
        return True
    
    def step5_extract_gene_sequences(self):
        """步骤5: 提取基因序列 (并行优化版本)"""
        self.log("🚀 开始步骤5: 提取基因序列")
        
        # 检查特异性基因列表
        specific_genes_file = self.genus_dir / "specific_genes.txt"
        if not specific_genes_file.exists():
            self.log("❌ 没有找到specific_genes.txt文件，请先运行步骤4", "ERROR")
            return False
        
        # 读取特异性基因列表
        with open(specific_genes_file, 'r') as f:
            specific_genes = [line.strip() for line in f if line.strip()]
        
        # 确定实际要处理的基因数量
        genes_to_process = specific_genes[:self.max_genes] if self.max_genes else specific_genes
        total_genes = len(specific_genes)
        processing_genes = len(genes_to_process)
        
        if self.max_genes and total_genes > self.max_genes:
            self.log(f"📋 总共有 {total_genes} 个特异性基因，限制处理前 {self.max_genes} 个")
        else:
            self.log(f"📋 处理所有 {processing_genes} 个特异性基因")
        
        # 创建序列提取目录
        extract_dir = self.genus_dir / "gene_sequences"
        extract_dir.mkdir(exist_ok=True)
        
        # 使用现有的extract脚本逻辑
        import pandas as pd
        from Bio import SeqIO
        
        # 查找gene_presence_absence.csv
        gene_presence_file = None
        for roary_dir in [self.roary_dir] + list(self.work_dir.glob("roary_*")):
            potential_file = roary_dir / "gene_presence_absence.csv"
            if potential_file.exists():
                gene_presence_file = potential_file
                break
        
        if not gene_presence_file:
            self.log("❌ 没有找到gene_presence_absence.csv文件", "ERROR")
            return False
        
        df = pd.read_csv(gene_presence_file, low_memory=False)
        
        # 获取样本列名
        sample_columns = [col for col in df.columns if col not in 
                         ['Gene', 'Non-unique Gene name', 'Annotation', 'No. isolates', 'No. sequences', 
                          'Avg sequences per isolate', 'Genome Fragment', 'Order within Fragment', 
                          'Accessory Fragment', 'Accessory Order with Fragment', 'QC', 'Min group size nuc', 
                          'Max group size nuc', 'Avg group size nuc']]
        
        # 查找ffn文件
        ffn_files = {}
        prokka_results = list(self.prokka_dir.glob("*/*.ffn"))
        
        for ffn_file in prokka_results:
            sample_name = ffn_file.parent.name
            # 尝试匹配样本名
            for col in sample_columns:
                if col in sample_name or sample_name in col:
                    ffn_files[col] = ffn_file
                    break
        
        self.log(f"📁 找到 {len(ffn_files)} 个ffn文件")
        
        if not ffn_files:
            self.log("❌ 没有找到对应的ffn文件", "ERROR")
            return False
        
        # 预读所有序列
        self.log("🔄 预加载序列文件...")
        seqdicts = {}
        
        def load_sequences(sample_ffn_pair):
            """并行加载序列文件"""
            sample, ffn_file = sample_ffn_pair
            try:
                seqdict = SeqIO.to_dict(SeqIO.parse(ffn_file, "fasta"))
                self.log(f"✅ 加载序列文件: {sample}")
                return sample, seqdict, True, None
            except Exception as e:
                self.log(f"❌ 加载序列文件失败 {sample}: {e}", "ERROR")
                return sample, None, False, str(e)
        
        # 并行加载序列文件
        max_workers = min(self.threads, len(ffn_files))
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_sample = {
                executor.submit(load_sequences, (sample, ffn_file)): sample 
                for sample, ffn_file in ffn_files.items()
            }
            
            for future in as_completed(future_to_sample):
                sample = future_to_sample[future]
                try:
                    sample_name, seqdict, success, error = future.result()
                    if success:
                        seqdicts[sample_name] = seqdict
                except Exception as exc:
                    self.log(f"❌ {sample} 序列加载出错: {exc}", "ERROR")
        
        # 并行提取特异性基因的序列
        def extract_gene_sequences(gene):
            """并行处理单个基因的序列提取"""
            try:
                gene_row = df[df['Gene'] == gene].iloc[0]
                seqs = []
                
                for sample in sample_columns:
                    if sample in seqdicts and pd.notna(gene_row[sample]):
                        gene_id = str(gene_row[sample]).strip()
                        if gene_id in seqdicts[sample]:
                            seq_record = seqdicts[sample][gene_id]
                            seq_record.id = f"{sample}_{gene_id}"
                            seq_record.description = f"{sample} {gene}"
                            seqs.append(seq_record)
                
                if len(seqs) >= 2:  # 至少需要2个序列才能进行比对 / At least 2 sequences needed for alignment
                    output_file = extract_dir / f"{gene}.fa"
                    SeqIO.write(seqs, output_file, "fasta")
                    self.log(f"✅ Extracted gene sequences: {gene} ({len(seqs)} sequences)")
                    return gene, len(seqs), True, None
                else:
                    self.log(f"⚠️  Insufficient gene sequences: {gene} (only {len(seqs)} sequences)")
                    return gene, len(seqs), False, "Insufficient sequence count"
            
            except Exception as e:
                self.log(f"❌ 提取基因序列失败 {gene}: {e}", "ERROR")
                return gene, 0, False, str(e)
        
        # 并行提取基因序列
        self.log("🔄 开始并行提取基因序列...")
        extracted_count = 0
        failed_genes = []
        
        max_workers = min(self.threads, processing_genes)
        self.log(f"🚀 使用 {max_workers} 个并行进程进行序列提取")
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # 提交所有序列提取任务
            future_to_gene = {
                executor.submit(extract_gene_sequences, gene): gene 
                for gene in genes_to_process
            }
            
            # 处理完成的任务
            for future in as_completed(future_to_gene):
                gene = future_to_gene[future]
                try:
                    gene_name, seq_count, success, error = future.result()
                    if success:
                        extracted_count += 1
                    else:
                        failed_genes.append(gene_name)
                except Exception as exc:
                    self.log(f"❌ {gene} 序列提取出错: {exc}", "ERROR")
                    failed_genes.append(gene)
        
        self.log(f"✅ 步骤5完成: 成功提取 {extracted_count}/{processing_genes} 个基因的序列")
        if failed_genes:
            self.log(f"❌ 失败的基因: {', '.join(failed_genes[:5])}{'...' if len(failed_genes) > 5 else ''}")
        
        return extracted_count > 0
    
    def step6_find_conserved_regions(self):
        """步骤6: 多序列比对并找到保守区域 (并行优化版本)"""
        self.log("🚀 开始步骤6: 多序列比对和保守区域识别")
        
        extract_dir = self.genus_dir / "gene_sequences"
        if not extract_dir.exists():
            self.log("❌ 没有找到基因序列目录，请先运行步骤5", "ERROR")
            return False
        
        # 创建比对结果目录
        alignment_dir = self.genus_dir / "alignments"
        alignment_dir.mkdir(exist_ok=True)
        
        # 获取所有基因序列文件
        gene_files = list(extract_dir.glob("*.fa"))
        if not gene_files:
            self.log("❌ 没有找到基因序列文件", "ERROR")
            return False
        
        self.log(f"📋 需要比对 {len(gene_files)} 个基因")
        
        # 检查mafft环境
        try:
            result = subprocess.run(["mafft", "--version"], capture_output=True, text=True)
            self.log("✅ MAFFT可用")
        except FileNotFoundError:
            self.log("❌ 未找到MAFFT程序，请安装MAFFT", "ERROR")
            return False
        
        # 并行执行多序列比对
        def align_single_gene(gene_file):
            """并行处理单个基因的比对"""
            gene_name = gene_file.stem
            alignment_file = alignment_dir / f"{gene_name}_aln.fa"
            
            try:
                self.log(f"🔄 正在比对: {gene_name}")
                
                # 运行MAFFT比对
                with open(gene_file, 'r') as input_f, open(alignment_file, 'w') as output_f:
                    result = subprocess.run(
                        ["mafft", "--auto", "--quiet", str(gene_file)],
                        stdout=output_f,
                        stderr=subprocess.PIPE,
                        text=True,
                        check=True
                    )
                
                self.log(f"✅ 完成比对: {gene_name}")
                return gene_name, alignment_file, True, None
                
            except Exception as e:
                self.log(f"❌ 比对失败 {gene_name}: {e}", "ERROR")
                return gene_name, None, False, str(e)
        
        # 计算合适的并行数量
        max_workers = min(self.threads, len(gene_files))
        self.log(f"🚀 使用 {max_workers} 个并行进程进行比对")
        
        # 执行并行比对
        aligned_count = 0
        alignment_files = []
        failed_genes = []
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # 提交所有比对任务
            future_to_gene = {
                executor.submit(align_single_gene, gene_file): gene_file 
                for gene_file in gene_files
            }
            
            # 处理完成的任务
            for future in as_completed(future_to_gene):
                gene_file = future_to_gene[future]
                try:
                    gene_name, alignment_file, success, error = future.result()
                    if success:
                        alignment_files.append(alignment_file)
                        aligned_count += 1
                    else:
                        failed_genes.append(gene_name)
                except Exception as exc:
                    gene_name = gene_file.stem
                    self.log(f"❌ {gene_name} 执行出错: {exc}", "ERROR")
                    failed_genes.append(gene_name)
        
        if aligned_count == 0:
            self.log("❌ 没有成功的比对结果", "ERROR")
            return False
        
        self.log(f"📊 成功比对 {aligned_count}/{len(gene_files)} 个基因")
        if failed_genes:
            self.log(f"❌ 失败的基因: {', '.join(failed_genes)}")
        
        # 第二步：并行查找保守区域
        self.log("🔍 开始查找保守区域...")
        
        def find_conserved_regions_for_gene(alignment_file):
            """并行处理单个基因的保守区域查找 - 寻找最长保守区域"""
            gene_name = alignment_file.stem.replace("_aln", "")
            conserved_regions = []
            
            try:
                # 查找保守区域 / Find conserved regions
                from Bio import AlignIO
                
                aln = AlignIO.read(alignment_file, "fasta")
                aln_len = aln.get_alignment_length()
                longest_conserved = None
                
                # 从最长到最短搜索保守区域，确保找到最长的 / Search from longest to shortest to ensure finding the longest
                # 搜索范围：从序列长度到50bp / Search range: from sequence length down to 50bp
                max_search_len = min(aln_len, 500)  # 最大搜索长度500bp / Maximum search length 500bp
                
                for win_len in range(max_search_len, 49, -1):  # 从长到短搜索 / Search from long to short
                    found_at_this_length = False
                    
                    for start in range(aln_len - win_len + 1):
                        window_seqs = [str(rec.seq[start:start+win_len]).upper() for rec in aln]
                        
                        # 检查是否完全一致且无gap / Check if completely identical and no gaps
                        if (len(set(window_seqs)) == 1 and  # 所有序列一致 / All sequences identical
                            "-" not in window_seqs[0]):    # 无gap / No gaps
                            
                            # 找到最长保守区域，立即返回 / Found longest conserved region, return immediately
                            longest_conserved = {
                                'gene': gene_name,
                                'position': f"{start+1}-{start+win_len}",
                                'length': win_len,
                                'sequence': window_seqs[0]
                            }
                            found_at_this_length = True
                            break
                    
                    # 如果在当前长度找到保守区域，这就是最长的，停止搜索 / If found at current length, this is the longest, stop searching
                    if found_at_this_length:
                        break
                
                if longest_conserved:
                    conserved_regions.append(longest_conserved)
                    length = longest_conserved['length']
                    if length >= 200:
                        self.log(f"🎯 发现长保守区域: {gene_name} ({length}bp)")
                    elif length >= 100:
                        self.log(f"✅ 发现保守区域: {gene_name} ({length}bp)")
                    elif length >= 80:
                        self.log(f"✅ 发现标准保守区域: {gene_name} ({length}bp)")
                    else:
                        self.log(f"⚠️  发现短保守区域: {gene_name} ({length}bp)")
                else:
                    conserved_regions.append({
                        'gene': gene_name,
                        'position': '无保守区域',
                        'length': 0,
                        'sequence': ''
                    })
                    self.log(f"❌ 无保守区域: {gene_name}")
                
                return gene_name, conserved_regions, longest_conserved is not None, None
                
            except Exception as e:
                self.log(f"❌ 分析保守区域失败 {gene_name}: {e}", "ERROR")
                return gene_name, [{
                    'gene': gene_name,
                    'position': '分析失败',
                    'length': 0,
                    'sequence': ''
                }], False, str(e)
        
        # 并行查找保守区域
        all_conserved_regions = []
        conserved_count = 0
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # 提交所有保守区域查找任务
            future_to_alignment = {
                executor.submit(find_conserved_regions_for_gene, alignment_file): alignment_file 
                for alignment_file in alignment_files
            }
            
            # 处理完成的任务
            for future in as_completed(future_to_alignment):
                alignment_file = future_to_alignment[future]
                try:
                    gene_name, regions, success, error = future.result()
                    all_conserved_regions.extend(regions)
                    if success and regions and regions[0]['length'] > 0:
                        conserved_count += 1
                except Exception as exc:
                    gene_name = alignment_file.stem.replace("_aln", "")
                    self.log(f"❌ {gene_name} 保守区域分析出错: {exc}", "ERROR")
        
        # 保存保守区域结果
        conserved_regions_file = self.genus_dir / "conserved_regions.txt"
        with open(conserved_regions_file, 'w') as conserved_f:
            conserved_f.write("Gene\tPosition\tLength\tSequence\n")
            for region in all_conserved_regions:
                conserved_f.write(f"{region['gene']}\t{region['position']}\t{region['length']}\t{region['sequence']}\n")
        
        self.log(f"✅ 步骤6完成: 找到 {conserved_count} 个保守区域")
        self.log(f"📁 比对结果保存到: {alignment_dir}")
        self.log(f"📁 保守区域结果保存到: {conserved_regions_file}")
        
        return conserved_count > 0
    
    def step7_primer3_design(self):
        """步骤7: 使用Primer3设计引物 (并行优化版本)"""
        self.log("🚀 开始步骤7: Primer3引物设计")
        
        conserved_regions_file = self.genus_dir / "conserved_regions.txt"
        if not conserved_regions_file.exists():
            self.log("❌ 没有找到保守区域文件，请先运行步骤6", "ERROR")
            return False
        
        # 检查primer3是否可用
        try:
            result = subprocess.run(["primer3_core", "--help"], capture_output=True, text=True)
            self.log("✅ Primer3可用")
        except FileNotFoundError:
            self.log("❌ 未找到primer3_core程序，请安装Primer3", "ERROR")
            return False
        
        # 创建primer3目录
        primer3_dir = self.genus_dir / "primer3_results"
        primer3_dir.mkdir(exist_ok=True)
        
        # 读取保守区域
        import pandas as pd
        conserved_df = pd.read_csv(conserved_regions_file, sep='\t')
        
        # 过滤有效的保守区域
        valid_regions = conserved_df[
            (conserved_df['Length'] >= 80) & 
            (conserved_df['Length'] <= 400) & 
            (conserved_df['Sequence'].notna())
        ]
        
        if len(valid_regions) == 0:
            self.log("❌ 没有找到有效的保守区域", "ERROR")
            return False
        
        self.log(f"📋 找到 {len(valid_regions)} 个有效保守区域")
        
        # 并行设计引物
        def design_primers_for_region(region_data):
            """并行处理单个保守区域的引物设计"""
            idx, row = region_data
            gene_name = row['Gene']
            sequence = row['Sequence']
            
            # 创建临时输入文件
            temp_input_file = primer3_dir / f"temp_input_{gene_name}_{idx}.txt"
            temp_output_file = primer3_dir / f"temp_output_{gene_name}_{idx}.txt"
            
            try:
                # 写入Primer3输入 (优化产物大小配置)
                with open(temp_input_file, 'w') as f:
                    f.write(f"SEQUENCE_ID={gene_name}_{idx}\n")
                    f.write(f"SEQUENCE_TEMPLATE={sequence}\n")
                    f.write("PRIMER_TASK=generic\n")
                    f.write("PRIMER_PICK_LEFT_PRIMER=1\n")
                    f.write("PRIMER_PICK_INTERNAL_OLIGO=0\n")
                    f.write("PRIMER_PICK_RIGHT_PRIMER=1\n")
                    f.write("PRIMER_OPT_SIZE=20\n")
                    f.write("PRIMER_MIN_SIZE=18\n")
                    f.write("PRIMER_MAX_SIZE=25\n")
                    f.write("PRIMER_MAX_NS_ACCEPTED=0\n")
                    # 优化产物大小配置 / Optimized product size configuration
                    f.write("PRIMER_PRODUCT_SIZE_RANGE=100-300 80-200 200-400\n")  # 多个范围选择 / Multiple range options
                    f.write("PRIMER_OPT_PRODUCT_SIZE=150\n")  # 理想产物大小 / Optimal product size
                    f.write("PRIMER_OPT_TM=60.0\n")
                    f.write("PRIMER_MIN_TM=55.0\n")
                    f.write("PRIMER_MAX_TM=65.0\n")
                    f.write("PRIMER_MIN_GC=30.0\n")
                    f.write("PRIMER_MAX_GC=70.0\n")
                    f.write("PRIMER_MAX_POLY_X=4\n")
                    # 移除有问题的配置路径
                    # f.write("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/usr/share/primer3/primer3_config/\n")
                    f.write("=\n")
                
                # 运行Primer3
                with open(temp_input_file, 'r') as input_f, open(temp_output_file, 'w') as output_f:
                    result = subprocess.run(
                        ["primer3_core"],
                        stdin=input_f,
                        stdout=output_f,
                        stderr=subprocess.PIPE,
                        text=True,
                        check=True
                    )
                
                # 读取输出结果
                with open(temp_output_file, 'r') as f:
                    output_content = f.read()
                
                # 清理临时文件
                temp_input_file.unlink(missing_ok=True)
                temp_output_file.unlink(missing_ok=True)
                
                self.log(f"✅ 完成引物设计: {gene_name}_{idx}")
                return gene_name, idx, output_content, True, None
                
            except subprocess.CalledProcessError as e:
                # 清理临时文件
                temp_input_file.unlink(missing_ok=True)
                temp_output_file.unlink(missing_ok=True)
                
                # 记录详细错误信息
                error_msg = f"Primer3 exit code: {e.returncode}"
                if e.stderr:
                    error_msg += f", stderr: {e.stderr[:200]}"
                
                self.log(f"❌ 引物设计失败 {gene_name}_{idx}: {error_msg}", "ERROR")
                return gene_name, idx, "", False, error_msg
                
            except Exception as e:
                # 清理临时文件
                temp_input_file.unlink(missing_ok=True)
                temp_output_file.unlink(missing_ok=True)
                
                self.log(f"❌ 引物设计失败 {gene_name}_{idx}: {e}", "ERROR")
                return gene_name, idx, "", False, str(e)
        
        # 准备并行任务数据
        region_tasks = [(idx, row) for idx, row in valid_regions.iterrows()]
        
        # 计算合适的并行数量
        max_workers = min(self.threads, len(region_tasks))
        self.log(f"🚀 使用 {max_workers} 个并行进程进行引物设计")
        
        # 执行并行引物设计
        all_outputs = []
        success_count = 0
        failed_regions = []
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # 提交所有引物设计任务
            future_to_region = {
                executor.submit(design_primers_for_region, region_data): region_data 
                for region_data in region_tasks
            }
            
            # 处理完成的任务
            for future in as_completed(future_to_region):
                region_data = future_to_region[future]
                try:
                    gene_name, idx, output_content, success, error = future.result()
                    if success:
                        all_outputs.append(output_content)
                        success_count += 1
                    else:
                        failed_regions.append(f"{gene_name}_{idx}")
                        # 记录第一个失败的详细错误
                        if len(failed_regions) == 1:
                            self.log(f"🔍 第一个失败的详细错误: {error}", "ERROR")
                except Exception as exc:
                    idx, row = region_data
                    gene_name = row['Gene']
                    self.log(f"❌ {gene_name}_{idx} 引物设计出错: {exc}", "ERROR")
                    failed_regions.append(f"{gene_name}_{idx}")
        
        # 合并所有输出结果
        primer3_output_file = primer3_dir / "primer3_output.txt"
        with open(primer3_output_file, 'w') as f:
            for output in all_outputs:
                f.write(output)
                if not output.endswith('\n'):
                    f.write('\n')
        
        self.log(f"📊 引物设计完成: 成功 {success_count}/{len(region_tasks)} 个区域")
        if failed_regions:
            self.log(f"❌ 失败的区域: {', '.join(failed_regions[:5])}{'...' if len(failed_regions) > 5 else ''}")
        
        if success_count > 0:
            self.log("✅ Primer3设计完成")
            self.log(f"📁 结果保存到: {primer3_output_file}")
            return True
        else:
            self.log("❌ 所有引物设计都失败了", "ERROR")
            return False
    
    def step8_parse_primer3_results(self):
        """Step 8: Parse Primer3 results using toolkit (并行优化版本)"""
        self.log("🚀 Starting Step 8: Parse Primer3 results")
        
        primer3_dir = self.genus_dir / "primer3_results"
        primer3_output_file = primer3_dir / "primer3_output.txt"
        
        if not primer3_output_file.exists():
            self.log("❌ Primer3 output file not found, please run step 7 first", "ERROR")
            return False
        
        # Use toolkit parser with parallel processing
        try:
            # 读取Primer3输出文件
            with open(primer3_output_file, 'r') as f:
                content = f.read()
            
            # 分割成独立的记录
            records = content.split('=\n')
            records = [record.strip() for record in records if record.strip()]
            
            self.log(f"📋 找到 {len(records)} 个Primer3记录需要解析")
            
            if not records:
                self.log("❌ 没有找到有效的Primer3记录", "ERROR")
                return False
            
            # 并行解析Primer3记录
            def parse_single_record(record):
                """并行处理单个Primer3记录的解析"""
                try:
                    lines = record.split('\n')
                    result = {}
                    
                    for line in lines:
                        if '=' in line:
                            key, value = line.split('=', 1)
                            result[key] = value
                    
                    # 提取基本信息
                    sequence_id = result.get('SEQUENCE_ID', '')
                    if not sequence_id:
                        return None, False, "缺少SEQUENCE_ID"
                    
                    # 解析引物信息
                    primers = []
                    
                    # 检查是否有引物返回
                    primer_returned = result.get('PRIMER_PAIR_NUM_RETURNED', '0')
                    if primer_returned == '0':
                        return {
                            'ID': sequence_id,
                            'Status': 'No primers found',
                            'Left_Primer': '',
                            'Right_Primer': '',
                            'Product_Size': 0,
                            'Left_Tm': 0,
                            'Right_Tm': 0,
                            'Left_GC': 0,
                            'Right_GC': 0
                        }, True, None
                    
                    # 解析第一个引物对
                    primer_data = {
                        'ID': sequence_id,
                        'Status': 'Success',
                        'Left_Primer': result.get('PRIMER_LEFT_0_SEQUENCE', ''),
                        'Right_Primer': result.get('PRIMER_RIGHT_0_SEQUENCE', ''),
                        'Product_Size': int(result.get('PRIMER_PAIR_0_PRODUCT_SIZE', 0)),
                        'Left_Tm': float(result.get('PRIMER_LEFT_0_TM', 0)),
                        'Right_Tm': float(result.get('PRIMER_RIGHT_0_TM', 0)),
                        'Left_GC': float(result.get('PRIMER_LEFT_0_GC_PERCENT', 0)),
                        'Right_GC': float(result.get('PRIMER_RIGHT_0_GC_PERCENT', 0))
                    }
                    
                    return primer_data, True, None
                    
                except Exception as e:
                    return None, False, str(e)
            
            # 并行解析所有记录
            max_workers = min(self.threads, len(records))
            self.log(f"🚀 使用 {max_workers} 个并行进程进行结果解析")
            
            parsed_results = []
            failed_count = 0
            
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                # 提交所有解析任务
                future_to_record = {
                    executor.submit(parse_single_record, record): i 
                    for i, record in enumerate(records)
                }
                
                # 处理完成的任务
                for future in as_completed(future_to_record):
                    record_idx = future_to_record[future]
                    try:
                        primer_data, success, error = future.result()
                        if success and primer_data:
                            parsed_results.append(primer_data)
                        else:
                            failed_count += 1
                            if error:
                                self.log(f"❌ 记录 {record_idx} 解析失败: {error}", "ERROR")
                    except Exception as exc:
                        failed_count += 1
                        self.log(f"❌ 记录 {record_idx} 解析出错: {exc}", "ERROR")
            
            # 保存解析结果
            parsed_results_file = primer3_dir / "parsed_primers.csv"
            
            if parsed_results:
                import pandas as pd
                df = pd.DataFrame(parsed_results)
                df.to_csv(parsed_results_file, index=False)
                
                self.log(f"📊 解析完成: 成功 {len(parsed_results)}/{len(records)} 个记录")
                if failed_count > 0:
                    self.log(f"❌ 失败记录数: {failed_count}")
                
                self.log("✅ Primer3 results parsing completed")
                self.log(f"📁 Parsed results saved to: {parsed_results_file}")
                return True
            else:
                self.log("❌ 没有成功解析任何记录", "ERROR")
                return False
            
        except Exception as e:
            self.log(f"❌ Failed to parse Primer3 results: {e}", "ERROR")
            return False
    
    def step9_rank_primers(self):
        """Step 9: Primer quality scoring and ranking (并行优化版本)"""
        self.log("🚀 Starting Step 9: Primer quality scoring and ranking")
        
        primer3_dir = self.genus_dir / "primer3_results"
        parsed_results_file = primer3_dir / "parsed_primers.csv"
        
        if not parsed_results_file.exists():
            self.log("❌ Parsed primer file not found, please run step 8 first", "ERROR")
            return False
        
        # 读取解析结果
        try:
            import pandas as pd
            df = pd.read_csv(parsed_results_file)
            
            if len(df) == 0:
                self.log("❌ 没有找到引物数据", "ERROR")
                return False
            
            self.log(f"📋 需要评分 {len(df)} 个引物")
            
            # 并行计算质量评分
            def calculate_quality_score(primer_data):
                """并行计算单个引物的质量评分"""
                try:
                    idx, row = primer_data
                    
                    # 检查引物是否成功设计
                    if row['Status'] != 'Success' or not row['Left_Primer'] or not row['Right_Primer']:
                        return idx, {
                            'Quality_Score': 0,
                            'Tm_Score': 0,
                            'GC_Score': 0,
                            'Length_Score': 0,
                            'Product_Score': 0,
                            'Overall_Grade': 'Failed'
                        }, True, None
                    
                    # 计算各项评分
                    scores = {}
                    
                    # 1. Tm评分 (理想范围: 55-65°C)
                    left_tm = float(row['Left_Tm'])
                    right_tm = float(row['Right_Tm'])
                    tm_diff = abs(left_tm - right_tm)
                    
                    # Tm一致性评分 (差异越小越好)
                    if tm_diff <= 2:
                        tm_consistency = 100
                    elif tm_diff <= 5:
                        tm_consistency = 80
                    elif tm_diff <= 10:
                        tm_consistency = 60
                    else:
                        tm_consistency = 40
                    
                    # Tm范围评分
                    avg_tm = (left_tm + right_tm) / 2
                    if 58 <= avg_tm <= 62:
                        tm_range = 100
                    elif 55 <= avg_tm <= 65:
                        tm_range = 80
                    elif 50 <= avg_tm <= 70:
                        tm_range = 60
                    else:
                        tm_range = 40
                    
                    scores['Tm_Score'] = (tm_consistency + tm_range) / 2
                    
                    # 2. GC含量评分 (理想范围: 40-60%)
                    left_gc = float(row['Left_GC'])
                    right_gc = float(row['Right_GC'])
                    avg_gc = (left_gc + right_gc) / 2
                    
                    if 45 <= avg_gc <= 55:
                        scores['GC_Score'] = 100
                    elif 40 <= avg_gc <= 60:
                        scores['GC_Score'] = 80
                    elif 35 <= avg_gc <= 65:
                        scores['GC_Score'] = 60
                    else:
                        scores['GC_Score'] = 40
                    
                    # 3. 引物长度评分 (理想长度: 18-25bp)
                    left_len = len(row['Left_Primer'])
                    right_len = len(row['Right_Primer'])
                    avg_len = (left_len + right_len) / 2
                    
                    if 19 <= avg_len <= 22:
                        scores['Length_Score'] = 100
                    elif 18 <= avg_len <= 25:
                        scores['Length_Score'] = 80
                    elif 16 <= avg_len <= 28:
                        scores['Length_Score'] = 60
                    else:
                        scores['Length_Score'] = 40
                    
                    # 4. 产物大小评分 (理想范围: 100-200bp)
                    product_size = int(row['Product_Size'])
                    if 120 <= product_size <= 180:
                        scores['Product_Score'] = 100
                    elif 100 <= product_size <= 200:
                        scores['Product_Score'] = 80
                    elif 80 <= product_size <= 250:
                        scores['Product_Score'] = 60
                    else:
                        scores['Product_Score'] = 40
                    
                    # 5. 计算总体质量评分 (加权平均)
                    weights = {
                        'Tm_Score': 0.3,
                        'GC_Score': 0.25,
                        'Length_Score': 0.2,
                        'Product_Score': 0.25
                    }
                    
                    quality_score = sum(scores[key] * weights[key] for key in weights)
                    scores['Quality_Score'] = round(quality_score, 2)
                    
                    # 6. 总体等级评定
                    if quality_score >= 85:
                        scores['Overall_Grade'] = 'Excellent'
                    elif quality_score >= 75:
                        scores['Overall_Grade'] = 'Good'
                    elif quality_score >= 65:
                        scores['Overall_Grade'] = 'Fair'
                    elif quality_score >= 50:
                        scores['Overall_Grade'] = 'Poor'
                    else:
                        scores['Overall_Grade'] = 'Very Poor'
                    
                    return idx, scores, True, None
                    
                except Exception as e:
                    return idx, {}, False, str(e)
            
            # 并行计算所有引物的质量评分
            max_workers = min(self.threads, len(df))
            self.log(f"🚀 使用 {max_workers} 个并行进程进行质量评分")
            
            all_scores = {}
            failed_count = 0
            
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                # 提交所有评分任务
                future_to_primer = {
                    executor.submit(calculate_quality_score, (idx, row)): idx 
                    for idx, row in df.iterrows()
                }
                
                # 处理完成的任务
                for future in as_completed(future_to_primer):
                    primer_idx = future_to_primer[future]
                    try:
                        idx, scores, success, error = future.result()
                        if success:
                            all_scores[idx] = scores
                        else:
                            failed_count += 1
                            if error:
                                self.log(f"❌ 引物 {idx} 评分失败: {error}", "ERROR")
                    except Exception as exc:
                        failed_count += 1
                        self.log(f"❌ 引物 {primer_idx} 评分出错: {exc}", "ERROR")
            
            # 将评分结果添加到DataFrame
            for score_key in ['Quality_Score', 'Tm_Score', 'GC_Score', 'Length_Score', 'Product_Score', 'Overall_Grade']:
                df[score_key] = df.index.map(lambda idx: all_scores.get(idx, {}).get(score_key, 0))
            
            # 按质量评分排序
            df_sorted = df.sort_values('Quality_Score', ascending=False)
            
            # 保存排序结果
            ranked_results_file = primer3_dir / "ranked_primers.csv"
            df_sorted.to_csv(ranked_results_file, index=False)
            
            self.log(f"📊 质量评分完成: 成功 {len(df) - failed_count}/{len(df)} 个引物")
            if failed_count > 0:
                self.log(f"❌ 失败引物数: {failed_count}")
            
            self.log("✅ Primer quality scoring and ranking completed")
            self.log(f"📁 Final results saved to: {ranked_results_file}")
            
            # Display top primers summary
            try:
                if len(df_sorted) > 0:
                    self.log(f"🏆 Total designed primers: {len(df_sorted)}")
                    self.log("🥇 Top 3 best primers:")
                    for i in range(min(3, len(df_sorted))):
                        row = df_sorted.iloc[i]
                        self.log(f"   {i+1}. {row['ID']} - Quality Score: {row.get('Quality_Score', 'N/A'):.2f} ({row.get('Overall_Grade', 'N/A')})")
                        
                    # Quality summary
                    grade_counts = df_sorted['Overall_Grade'].value_counts()
                    self.log(f"📊 Quality Summary:")
                    for grade, count in grade_counts.items():
                        self.log(f"   {grade}: {count}")
                        
            except Exception as e:
                self.log(f"⚠️  Error displaying results summary: {e}")
            
            return True
            
        except Exception as e:
            self.log(f"❌ Failed to rank primers: {e}", "ERROR")
            return False
    
    def expand_gene_selection(self):
        """在引物设计失败时扩展基因选择范围"""
        self.log("🔄 引物设计效果不佳，尝试扩展基因选择范围")
        
        # 检查是否有保存的所有特异性基因信息
        detailed_genes_file = self.genus_dir / "specific_genes_detailed.csv"
        if not detailed_genes_file.exists():
            self.log("❌ 没有找到详细基因信息文件", "ERROR")
            return False
        
        # 检查roary结果
        gene_presence_file = None
        for roary_dir in [self.roary_dir] + list(self.work_dir.glob("roary_*")):
            potential_file = roary_dir / "gene_presence_absence.csv"
            if potential_file.exists():
                gene_presence_file = potential_file
                break
        
        if not gene_presence_file:
            self.log("❌ 没有找到gene_presence_absence.csv文件", "ERROR")
            return False
        
        import pandas as pd
        df = pd.read_csv(gene_presence_file, low_memory=False)
        
        # 获取样本列名
        sample_columns = [col for col in df.columns if col not in 
                         ['Gene', 'Non-unique Gene name', 'Annotation', 'No. isolates', 'No. sequences', 
                          'Avg sequences per isolate', 'Genome Fragment', 'Order within Fragment', 
                          'Accessory Fragment', 'Accessory Order with Fragment', 'QC', 'Min group size nuc', 
                          'Max group size nuc', 'Avg group size nuc']]
        
        # 重新分类样本
        target_prefix = self.genus[:3].capitalize() + "_"
        target_samples = [sample for sample in sample_columns if sample.startswith(target_prefix)]
        outgroup_samples = [sample for sample in sample_columns if not sample.startswith(target_prefix)]
        
        # 重新筛选所有特异性基因
        all_specific_genes = []
        for idx, row in df.iterrows():
            gene_name = row['Gene']
            
            target_present = all(pd.notna(row[sample]) and str(row[sample]).strip() != '' 
                               for sample in target_samples)
            
            if outgroup_samples:
                outgroup_absent = all(pd.isna(row[sample]) or str(row[sample]).strip() == '' 
                                    for sample in outgroup_samples)
            else:
                outgroup_absent = True
            
            if target_present and outgroup_absent:
                all_specific_genes.append(gene_name)
        
        # 更新基因列表
        specific_genes_file = self.genus_dir / "specific_genes.txt"
        with open(specific_genes_file, 'w') as f:
            for gene in all_specific_genes:
                f.write(f"{gene}\n")
        
        # 更新详细信息
        selected_df = df[df['Gene'].isin(all_specific_genes)][['Gene', 'Annotation']].copy()
        selected_df.to_csv(detailed_genes_file, index=False)
        
        self.log(f"✅ 已扩展到所有 {len(all_specific_genes)} 个特异性基因")
        self.log(f"📁 更新的基因列表保存到: {specific_genes_file}")
        
        return True
    
    def run_pipeline(self, start_step=1):
        """运行完整流程"""
        self.log(f"🎯 开始运行基因组分析流程 - {self.genus}")
        self.log(f"📍 工作目录: {self.work_dir}")
        self.log(f"🔢 从步骤 {start_step} 开始")
        
        # 创建工作目录
        self.genus_dir.mkdir(exist_ok=True)
        
        steps = {
            1: ("下载基因组", self.step1_download),
            2: ("Prokka注释", self.step2_prokka),
            3: ("Roary泛基因组分析", self.step3_roary),
            4: ("筛选目标特异性基因", self.step4_find_specific_genes),
            5: ("提取基因序列并进行多序列比对", self.step5_extract_gene_sequences),
            6: ("多序列比对并找到保守区域", self.step6_find_conserved_regions),
            7: ("使用Primer3设计引物", self.step7_primer3_design),
            8: ("解析Primer3结果", self.step8_parse_primer3_results),
            9: ("引物质量评分和排序", self.step9_rank_primers)
        }
        
        for step_num in range(start_step, 10):
            step_name, step_func = steps[step_num]
            self.log(f"\n{'='*50}")
            self.log(f"执行步骤 {step_num}: {step_name}")
            self.log(f"{'='*50}")
            
            if not step_func():
                self.log(f"❌ 步骤 {step_num} 失败，流程终止", "ERROR")
                return False
        
        self.log("\n🎉 完整流程执行成功！")
        self.log("📁 结果文件位置:")
        self.log(f"   - 基因组文件: {self.data_dir}")
        
        # 显示基因组summary文件
        summary_files = []
        target_summary = self.genus_dir / f"{self.genus_safe}_genome_summary.csv"
        if target_summary.exists():
            summary_files.append(f"   - 目标基因组summary: {target_summary}")
        
        for outgroup, outgroup_safe in zip(self.outgroup_genera, self.outgroup_genera_safe):
            outgroup_dir = self.work_dir / outgroup_safe
            outgroup_summary = outgroup_dir / f"{outgroup_safe}_genome_summary.csv"
            if outgroup_summary.exists():
                summary_files.append(f"   - 外群基因组summary: {outgroup_summary}")
        
        if summary_files:
            for summary_file in summary_files:
                self.log(summary_file)
        
        self.log(f"   - Prokka注释: {self.prokka_dir}")
        
        # 检查roary结果目录
        roary_result_dir = None
        for roary_dir in [self.roary_dir] + list(self.work_dir.glob("roary_*")):
            if (roary_dir / "gene_presence_absence.csv").exists():
                roary_result_dir = roary_dir
                break
        if roary_result_dir:
            self.log(f"   - Roary结果: {roary_result_dir}")
        
        # 检查引物设计结果
        if (self.genus_dir / "specific_genes.txt").exists():
            self.log(f"   - 特异性基因: {self.genus_dir / 'specific_genes.txt'}")
        if (self.genus_dir / "gene_sequences").exists():
            self.log(f"   - 基因序列: {self.genus_dir / 'gene_sequences'}")
        if (self.genus_dir / "conserved_regions.txt").exists():
            self.log(f"   - 保守区域: {self.genus_dir / 'conserved_regions.txt'}")
        if (self.genus_dir / "primer3_results" / "ranked_primers.csv").exists():
            self.log(f"   - 最终引物: {self.genus_dir / 'primer3_results' / 'ranked_primers.csv'}")
        
        return True

def main():
    parser = argparse.ArgumentParser(
        description="完整的基因组分析流程：下载 -> Prokka -> Roary -> 引物设计",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法:
  # 完整流程 (下载 -> Prokka -> Roary -> 引物设计，默认快速模式)
  python complete_pipeline.py Intestinibacter --outgroup Clostridium Bacteroides
  
  # 使用高质量模式 (适合最终分析和引物设计)
  python complete_pipeline.py Intestinibacter --high-quality
  
  # 从步骤2开始（跳过下载）
  python complete_pipeline.py Intestinibacter --start-step 2
  
  # 只运行Roary分析 (快速模式)
  python complete_pipeline.py Intestinibacter --start-step 3 --threads 50
  
  # 从引物设计开始 (假设已有Roary结果)
  python complete_pipeline.py Intestinibacter --start-step 4
  
  # 只运行引物质量评分 (假设已有解析结果)
  python complete_pipeline.py Intestinibacter --start-step 9

基因选择策略:
  步骤4会智能选择特异性基因用于引物设计：
  1. 优先选择hypothetical protein基因(≥5个)
  2. 如果hypothetical protein不足，则使用所有特异性基因
  3. 如果引物设计效果不佳，可以手动扩展基因范围：
     - 修改specific_genes.txt文件
     - 或者重新运行步骤5-9

外群识别:
  自动识别：目标genus前缀开头的为目标基因组，其他为外群
  例如：Terrisporobacter -> Ter_开头为目标，其他(如Int_)为外群
  
步骤说明:
  1: 下载基因组
  2: Prokka注释
  3: Roary泛基因组分析
  4: 筛选目标特异性基因 (优先hypothetical protein)
  5: 提取基因序列
  6: 多序列比对和保守区域识别
  7: Primer3引物设计
  8: 解析Primer3结果
  9: 引物质量评分和排序
        """
    )
    
    parser.add_argument("genus", help="目标genus名称")
    parser.add_argument("--outgroup", nargs="*", help="外群genus名称列表")
    parser.add_argument("--level", default="Complete Genome", help="组装级别 (默认: 'Complete Genome')")
    parser.add_argument("--threads", type=int, default=4, help="并行线程数 (默认: 4)")
    parser.add_argument("--start-step", type=int, choices=[1, 2, 3, 4, 5, 6, 7, 8, 9], default=1,
                       help="开始步骤: 1=下载, 2=Prokka, 3=Roary, 4=筛选目标特异性基因, 5=提取基因序列并进行多序列比对, 6=多序列比对并找到保守区域, 7=使用Primer3设计引物, 8=解析Primer3结果, 9=引物质量评分和排序 (默认: 1)")
    parser.add_argument("--fast", action="store_true", default=True,
                       help="使用快速模式 (-n) 进行Roary分析 (默认: True)")
    parser.add_argument("--high-quality", action="store_true", 
                       help="使用高质量模式 (--mafft) 进行Roary分析，覆盖 --fast")
    
    args = parser.parse_args()
    
    # 确定使用哪种模式
    fast_mode = not args.high_quality  # 如果指定了high_quality，就不使用fast模式
    
    # 创建流程实例
    pipeline = GenomePipeline(
        genus=args.genus,
        outgroup_genera=args.outgroup,
        level=args.level,
        threads=args.threads,
        fast_mode=fast_mode
    )
    
    # 运行流程
    success = pipeline.run_pipeline(start_step=args.start_step)
    
    if not success:
        sys.exit(1)

if __name__ == "__main__":
    main() 