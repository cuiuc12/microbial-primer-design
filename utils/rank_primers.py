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
        print("📊 引物质量评分工具 / Primer Quality Ranking Tool")
        print(f"📥 输入文件 / Input file: {args.input_file}")
        print(f"📤 输出文件 / Output file: {args.output_file}")
        
        # 检查输入文件
        if not Path(args.input_file).exists():
            print(f"❌ 输入文件不存在 / Input file not found: {args.input_file}")
            sys.exit(1)
        
        # 质量评分和排序
        ranker = PrimerQualityRanker()
        df = ranker.rank_primers(args.input_file, args.output_file)
        
        print(f"✅ 成功评分 {len(df)} 对引物 / Successfully ranked {len(df)} primer pairs")
        print(f"📁 结果已保存到 / Results saved to: {args.output_file}")
        
        if len(df) > 0:
            print(f"\n🏆 最佳引物 / Top primer: {df.iloc[0]['ID']} (分数/Score: {df.iloc[0]['Quality_Score']:.2f})")
        
        if args.summary or args.verbose:
            summary = ranker.get_quality_summary(df)
            print("\n📈 质量摘要 / Quality Summary:")
            print(f"   总引物数 / Total primers: {summary.get('total_primers', 0)}")
            print(f"   高质量 (≥80) / High quality: {summary.get('high_quality', 0)}")
            print(f"   中等质量 (60-79) / Medium quality: {summary.get('medium_quality', 0)}")
            print(f"   低质量 (<60) / Low quality: {summary.get('low_quality', 0)}")
            print(f"   平均分数 / Average score: {summary.get('avg_quality_score', 0):.2f}")
            
            if args.verbose:
                print(f"   平均Tm(左) / Avg Tm(left): {summary.get('avg_tm_left', 0):.1f}°C")
                print(f"   平均Tm(右) / Avg Tm(right): {summary.get('avg_tm_right', 0):.1f}°C")
                print(f"   平均GC(左) / Avg GC(left): {summary.get('avg_gc_left', 0):.1f}%")
                print(f"   平均GC(右) / Avg GC(right): {summary.get('avg_gc_right', 0):.1f}%")
                print(f"   平均产物大小 / Avg product size: {summary.get('avg_product_size', 0):.1f}bp")
        
    except Exception as e:
        print(f"❌ 评分失败 / Ranking failed: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main() 