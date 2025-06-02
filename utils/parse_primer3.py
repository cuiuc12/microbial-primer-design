#!/usr/bin/env python3
"""
Enhanced Primer3结果解析工具 / Enhanced Primer3 Result Parser Tool

完整的Primer3输出解析工具，支持二聚体分析 / Complete Primer3 output parser with dimer analysis support
"""

import sys
import argparse
from pathlib import Path

# Add toolkit to path
toolkit_path = Path(__file__).parent.parent / "primer_design_toolkit"
sys.path.insert(0, str(toolkit_path))

from primer3_parser import Primer3Parser


def main():
    """主函数 / Main function"""
    
    parser = argparse.ArgumentParser(
        description="增强版Primer3输出文件解析器 / Enhanced Primer3 output file parser",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例 / Usage Examples:
  python parse_primer3.py primer3_output.txt parsed_primers.csv
  python parse_primer3.py primer3_output.txt parsed_primers.csv --verbose
  
功能说明 / Function Description:
  将Primer3的原始输出转换为结构化的CSV格式，包括：
  Convert Primer3 raw output to structured CSV format, including:
  
  ✅ 基本引物参数 / Basic primer parameters:
     - 引物序列、Tm值、GC含量、长度 / Primer sequences, Tm, GC content, length
     - 产物大小和位置信息 / Product size and position information
     
  ✅ 二聚体分析参数 / Dimer analysis parameters:
     - 自身二聚体分析 (Self_Any, Self_End) / Self-dimer analysis
     - 发夹结构分析 (Hairpin) / Hairpin structure analysis
     - 引物对互补性分析 (Compl_Any, Compl_End) / Primer pair complementarity
     
  ✅ 衍生质量指标 / Derived quality metrics:
     - Tm差异、GC差异、长度差异 / Tm diff, GC diff, length diff
        """
    )
    
    parser.add_argument("input_file", help="Primer3输出文件 / Primer3 output file")
    parser.add_argument("output_file", help="输出CSV文件 / Output CSV file")
    parser.add_argument("--verbose", "-v", action="store_true", 
                       help="详细输出，显示解析统计 / Verbose output with parsing statistics")
    parser.add_argument("--summary", "-s", action="store_true",
                       help="显示详细摘要信息 / Show detailed summary information")
    
    args = parser.parse_args()
    
    try:
        print("🔧 Enhanced Primer3 Result Parser")
        print("=" * 50)
        print(f"📥 Input file: {args.input_file}")
        print(f"📤 Output file: {args.output_file}")
        print("")
        
        # 检查输入文件 / Check input file
        if not Path(args.input_file).exists():
            print(f"❌ Input file not found: {args.input_file}")
            sys.exit(1)
        
        # 解析文件 / Parse file
        print("🔄 Parsing Primer3 output...")
        parser_obj = Primer3Parser()
        df = parser_obj.parse_file(args.input_file, args.output_file)
        
        print(f"\n✅ Successfully parsed {len(df)} primer pairs")
        print(f"📁 Results saved to: {args.output_file}")
        
        if len(df) > 0 and (args.verbose or args.summary):
            print("\n" + "=" * 50)
            print("📊 Parsing Summary:")
            print("=" * 50)
            
            # 基本统计 / Basic statistics
            print(f"📋 Basic Statistics:")
            print(f"   Total primer pairs: {len(df)}")
            print(f"   Unique genes: {df['Gene'].nunique()}")
            print(f"   Average product size: {df['Product_Size'].mean():.1f} bp")
            print(f"   Product size range: {df['Product_Size'].min()}-{df['Product_Size'].max()} bp")
            
            # 引物参数统计 / Primer parameter statistics
            print(f"\n🧬 Primer Parameters:")
            print(f"   Average left Tm: {df['Left_Tm'].mean():.1f}°C")
            print(f"   Average right Tm: {df['Right_Tm'].mean():.1f}°C")
            print(f"   Average left GC: {df['Left_GC%'].mean():.1f}%")
            print(f"   Average right GC: {df['Right_GC%'].mean():.1f}%")
            
            # 二聚体分析统计 / Dimer analysis statistics
            if 'Left_Self_Any' in df.columns:
                print(f"\n🔗 Dimer Analysis:")
                print(f"   Average left self-dimer: {df['Left_Self_Any'].mean():.2f}")
                print(f"   Average right self-dimer: {df['Right_Self_Any'].mean():.2f}")
                if 'Pair_Compl_Any' in df.columns:
                    print(f"   Average pair complementarity: {df['Pair_Compl_Any'].mean():.2f}")
                
                # 质量分级 / Quality grading
                good_dimers = len(df[(df['Left_Self_Any'] <= 8) & (df['Right_Self_Any'] <= 8)])
                print(f"   Good dimer performance (≤8): {good_dimers}/{len(df)} ({good_dimers/len(df)*100:.1f}%)")
            
            # 列信息 / Column information
            print(f"\n📝 Output Columns:")
            basic_cols = [c for c in df.columns if not any(x in c for x in ['Self', 'Hairpin', 'Compl'])]
            dimer_cols = [c for c in df.columns if any(x in c for x in ['Self', 'Hairpin', 'Compl'])]
            
            print(f"   Basic parameters: {len(basic_cols)} columns")
            if args.verbose:
                print(f"      {', '.join(basic_cols[:5])}{'...' if len(basic_cols) > 5 else ''}")
            
            print(f"   Dimer parameters: {len(dimer_cols)} columns")
            if args.verbose and dimer_cols:
                print(f"      {', '.join(dimer_cols[:5])}{'...' if len(dimer_cols) > 5 else ''}")
            
            # 显示示例结果 / Show sample results
            if args.summary and len(df) > 0:
                print(f"\n🔍 Sample Results (Top 3):")
                sample_cols = ['ID', 'Product_Size', 'Left_Tm', 'Right_Tm', 'Left_Self_Any', 'Right_Self_Any']
                available_cols = [c for c in sample_cols if c in df.columns]
                if available_cols:
                    print(df[available_cols].head(3).to_string(index=False))
        
        print(f"\n🎉 Parsing completed successfully!")
        
    except Exception as e:
        print(f"❌ Parsing failed: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main() 