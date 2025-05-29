#!/usr/bin/env python3

import sys
import argparse
from pathlib import Path

sys.path.append(str(Path(__file__).parent))

from primer_design_toolkit.primer_pipeline import GenomePipeline

def main():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument("genus", help="目标genus名称 / Target genus name")
    parser.add_argument("--outgroup", nargs="*", help="外群genus名称列表 / Outgroup genus name list")
    parser.add_argument("--level", default="Complete Genome", help="组装级别 / Assembly level (default: 'Complete Genome')")
    parser.add_argument("--threads", type=int, default=4, 
                       help="并行线程数 / Number of parallel threads (default: 4, recommended: 20-50)")
    parser.add_argument("--start-step", type=int, choices=[1, 2, 3, 4, 5, 6, 7, 8, 9], default=1,
                       help="开始步骤 / Starting step (supports resume): 1=download, 2=Prokka, 3=Roary, 4=filter genes, 5=extract sequences, 6=alignment, 7=primer design, 8=parse results, 9=quality scoring (default: 1)")
    parser.add_argument("--fast", action="store_true", default=True,
                       help="使用快速模式 / Use fast mode (-n) for Roary analysis (default: True)")
    parser.add_argument("--high-quality", action="store_true", 
                       help="使用高质量模式 / Use high quality mode (--mafft) for Roary analysis, overrides --fast")
    parser.add_argument("--max-genes", type=int, default=None,
                       help="最大处理基因数量 / Maximum number of genes to process (default: process all genes, recommended: 20-100)")
    parser.add_argument("--assembly-summary", type=str, default=None,
                       help="assembly_summary.txt文件路径 / Path to assembly_summary.txt file (default: auto-detect)")
    
    args = parser.parse_args()
    
    # 确定使用哪种模式 / Determine which mode to use
    fast_mode = not args.high_quality  # 如果指定了high_quality，就不使用fast模式 / If high_quality is specified, don't use fast mode
    
    # 显示配置信息 / Display configuration information
    print("🚀 Microbial Primer Design Pipeline (Fully Parallelized Optimized Version)")
    print("=" * 70)
    print(f"🎯 Target genus: {args.genus}")
    if args.outgroup:
        print(f"🔗 Outgroups: {', '.join(args.outgroup)}")
    print(f"🧵 Parallel threads: {args.threads}")
    print(f"📊 Max genes: {args.max_genes or 'All'}")
    print(f"⚡ Analysis mode: {'High quality mode (--mafft)' if not fast_mode else 'Fast mode (-n)'}")
    print(f"🔄 Starting step: {args.start_step}")
    print("=" * 70)
    
    # 创建流程实例 / Create pipeline instance
    pipeline = GenomePipeline(
        genus=args.genus,
        outgroup_genera=args.outgroup,
        level=args.level,
        threads=args.threads,
        fast_mode=fast_mode,
        assembly_summary_path=args.assembly_summary,
        max_genes=args.max_genes
    )
    
    # 运行流程 / Run pipeline
    success = pipeline.run_pipeline(start_step=args.start_step)
    
    if success:
        print("\n🎉 Pipeline executed successfully!")
    else:
        print("\n❌ Pipeline execution failed")
        sys.exit(1)

if __name__ == "__main__":
    main() 
