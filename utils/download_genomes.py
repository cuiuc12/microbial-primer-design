#!/usr/bin/env python3
"""
基因组下载工具 / Genome Download Tool
"""

import sys
import argparse
from pathlib import Path

# Add toolkit to path
toolkit_path = Path(__file__).parent.parent / "primer_design_toolkit"
sys.path.insert(0, str(toolkit_path))

from genome_downloader import GenomeDownloader


def main():
    """主函数 / Main function"""
    
    parser = argparse.ArgumentParser(
        description="从NCBI下载细菌基因组 / Download bacterial genomes from NCBI",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例 / Usage Examples:
  # 下载单个属 / Download single genus
  python download_genomes.py Terrisporobacter
  
  # 下载目标属和外群 / Download target genus with outgroups
  python download_genomes.py Terrisporobacter --outgroup Intestinibacter Clostridium
  
  # 自定义组装级别 / Custom assembly level
  python download_genomes.py Terrisporobacter --level "Complete Genome"
  
功能说明 / Function Description:
  使用NCBI datasets工具下载并整理细菌基因组文件
  Download and organize bacterial genome files using NCBI datasets tool
        """
    )
    
    parser.add_argument("genus", help="目标属名 / Target genus name")
    parser.add_argument("--outgroup", nargs="*", 
                       help="外群属名列表 / Outgroup genus names")
    parser.add_argument("--level", default="Complete Genome", 
                       help="组装级别 / Assembly level (default: 'Complete Genome')")
    parser.add_argument("--threads", type=int, default=4, 
                       help="线程数 / Number of threads (default: 4)")
    parser.add_argument("--work-dir", 
                       help="工作目录 / Working directory (default: current)")
    parser.add_argument("--verbose", "-v", action="store_true",
                       help="详细输出 / Verbose output")
    
    args = parser.parse_args()
    
    try:
        print("📥 基因组下载工具 / Genome Download Tool")
        print(f"🎯 目标属 / Target genus: {args.genus}")
        if args.outgroup:
            print(f"🔗 外群 / Outgroups: {', '.join(args.outgroup)}")
        print(f"📊 组装级别 / Assembly level: {args.level}")
        print(f"⚙️  线程数 / Threads: {args.threads}")
        
        # 创建下载器
        downloader = GenomeDownloader(args.work_dir)
        
        # 下载基因组
        if args.outgroup:
            print("\n🚀 开始下载目标属和外群基因组 / Starting download of target and outgroup genomes")
            success = downloader.download_with_outgroup(
                args.genus, args.outgroup, args.level, args.threads
            )
        else:
            print(f"\n🚀 开始下载 {args.genus} 基因组 / Starting download of {args.genus} genomes")
            success = downloader.download_genus(
                args.genus, args.level, args.threads
            )
        
        if success:
            print("\n✅ 下载完成！/ Download completed successfully!")
            
            # 显示下载结果
            work_dir = Path(args.work_dir) if args.work_dir else Path.cwd()
            genus_dir = work_dir / args.genus
            target_dir = genus_dir / "data" / args.genus[:3].capitalize()
            outgroup_dir = genus_dir / "data" / "outgroup"
            
            if target_dir.exists():
                target_files = list(target_dir.glob("*.fna"))
                print(f"📁 目标基因组 / Target genomes: {len(target_files)} 个文件 / files")
                if args.verbose:
                    for f in target_files[:3]:  # 显示前3个
                        print(f"   - {f.name}")
                    if len(target_files) > 3:
                        print(f"   - ... 还有 {len(target_files)-3} 个文件 / and {len(target_files)-3} more")
            
            if outgroup_dir.exists():
                outgroup_files = list(outgroup_dir.glob("*.fna"))
                if outgroup_files:
                    print(f"📁 外群基因组 / Outgroup genomes: {len(outgroup_files)} 个文件 / files")
                    if args.verbose:
                        for f in outgroup_files:
                            print(f"   - {f.name}")
            
            print(f"\n📂 结果目录 / Result directory: {genus_dir}")
            
        else:
            print("❌ 下载失败 / Download failed")
            sys.exit(1)
            
    except KeyboardInterrupt:
        print("\n⚠️  用户中断下载 / User interrupted download")
        sys.exit(1)
    except Exception as e:
        print(f"❌ 下载过程中发生错误 / Error during download: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main() 