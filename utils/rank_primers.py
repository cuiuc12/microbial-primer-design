#!/usr/bin/env python3
"""
引物质量评分工具 / Primer Quality Ranking Tool
"""

import sys
import argparse
from pathlib import Path

# Add toolkit to path
toolkit_path = Path(__file__).parent.parent / "primer_design_toolkit"
sys.path.insert(0, str(toolkit_path))

from quality_ranker import PrimerQualityRanker


def main():
    """主函数 / Main function"""
    
    parser = argparse.ArgumentParser(
        description="引物质量评分和排序 / Primer quality scoring and ranking",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例 / Usage Examples:
  python rank_primers.py parsed_primers.csv ranked_primers.csv
  python rank_primers.py parsed_primers.csv ranked_primers.csv --summary
  
评分标准 / Scoring Criteria:
  - 熔解温度 (Tm): 最适60°C / Melting temperature: optimal 60°C
  - GC含量: 最适50% / GC content: optimal 50%  
  - 产物大小: 最适100-200bp / Product size: optimal 100-200bp
  - Tm差异: 越小越好 / Tm difference: smaller is better
  - 引物长度: 最适18-25nt / Primer length: optimal 18-25nt
        """
    )
    
    parser.add_argument("input_file", help="解析后的引物CSV文件 / Parsed primer CSV file")
    parser.add_argument("output_file", help="排序后的输出文件 / Ranked output file")
    parser.add_argument("--summary", action="store_true", 
                       help="显示质量摘要 / Show quality summary")
    parser.add_argument("--verbose", "-v", action="store_true",
                       help="详细输出 / Verbose output")
    
    args = parser.parse_args()
    
    try:
        print("📊 Primer Quality Ranking Tool")
        print(f"📥 Input file: {args.input_file}")
        print(f"📤 Output file: {args.output_file}")
        
        # 检查输入文件 / Check input file
        if not Path(args.input_file).exists():
            print(f"❌ Input file not found: {args.input_file}")
            sys.exit(1)
        
        # 质量评分和排序 / Quality scoring and ranking
        ranker = PrimerQualityRanker()
        df = ranker.rank_primers(args.input_file, args.output_file)
        
        print(f"✅ Successfully ranked {len(df)} primer pairs")
        print(f"📁 Results saved to: {args.output_file}")
        
        if len(df) > 0:
            print(f"\n🏆 Top primer: {df.iloc[0]['ID']} (Score: {df.iloc[0]['Quality_Score']:.2f})")
        
        if args.summary or args.verbose:
            summary = ranker.get_quality_summary(df)
            print("\n📈 Quality Summary:")
            print(f"   Total primers: {summary.get('total_primers', 0)}")
            print(f"   High quality (≥80): {summary.get('high_quality', 0)}")
            print(f"   Medium quality (60-79): {summary.get('medium_quality', 0)}")
            print(f"   Low quality (<60): {summary.get('low_quality', 0)}")
            print(f"   Average score: {summary.get('avg_quality_score', 0):.2f}")
            
            if args.verbose:
                print(f"   Avg Tm(left): {summary.get('avg_tm_left', 0):.1f}°C")
                print(f"   Avg Tm(right): {summary.get('avg_tm_right', 0):.1f}°C")
                print(f"   Avg GC(left): {summary.get('avg_gc_left', 0):.1f}%")
                print(f"   Avg GC(right): {summary.get('avg_gc_right', 0):.1f}%")
                print(f"   Avg product size: {summary.get('avg_product_size', 0):.1f}bp")
        
    except Exception as e:
        print(f"❌ Ranking failed: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main() 