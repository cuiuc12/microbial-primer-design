#!/usr/bin/env python3
"""
Primer3结果解析工具 / Primer3 Result Parser Tool
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
        description="解析Primer3输出文件 / Parse Primer3 output files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例 / Usage Examples:
  python parse_primer3.py primer3_output.txt parsed_primers.csv
  
功能说明 / Function Description:
  将Primer3的原始输出转换为结构化的CSV格式
  Convert Primer3 raw output to structured CSV format
        """
    )
    
    parser.add_argument("input_file", help="Primer3输出文件 / Primer3 output file")
    parser.add_argument("output_file", help="输出CSV文件 / Output CSV file")
    parser.add_argument("--verbose", "-v", action="store_true", 
                       help="详细输出 / Verbose output")
    
    args = parser.parse_args()
    
    try:
        print("🔧 Primer3结果解析工具 / Primer3 Result Parser")
        print(f"📥 输入文件 / Input file: {args.input_file}")
        print(f"📤 输出文件 / Output file: {args.output_file}")
        
        # 检查输入文件
        if not Path(args.input_file).exists():
            print(f"❌ 输入文件不存在 / Input file not found: {args.input_file}")
            sys.exit(1)
        
        # 解析文件
        parser_obj = Primer3Parser()
        df = parser_obj.parse_file(args.input_file, args.output_file)
        
        print(f"✅ 成功解析 {len(df)} 对引物 / Successfully parsed {len(df)} primer pairs")
        print(f"📁 结果已保存到 / Results saved to: {args.output_file}")
        
        if args.verbose and len(df) > 0:
            print("\n📊 解析结果摘要 / Parsing Summary:")
            print(f"   引物对数量 / Number of primer pairs: {len(df)}")
            print(f"   基因数量 / Number of genes: {df['Gene'].nunique()}")
            print(f"   平均产物大小 / Average product size: {df['Product_Size'].mean():.1f} bp")
        
    except Exception as e:
        print(f"❌ 解析失败 / Parsing failed: {e}")
        if args.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main() 