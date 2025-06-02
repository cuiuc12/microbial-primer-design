#!/usr/bin/env python3
"""
Enhanced Primer3 Result Parser / 增强版 Primer3 结果解析器
Parse Primer3 output files and extract comprehensive primer information including dimer analysis
解析Primer3输出文件并提取包括二聚体分析在内的全面引物信息
"""

import pandas as pd
import re
from pathlib import Path
import logging


class Primer3Parser:
    """
    Enhanced Parser for Primer3 output files / 增强版 Primer3 输出文件解析器
    Extracts primer pair information including dimer parameters and converts to structured format
    提取包括二聚体参数在内的引物对信息并转换为结构化格式
    """
    
    def __init__(self):
        self.primers = []
        self.logger = logging.getLogger(__name__)
        
    def parse_file(self, input_file, output_file=None):
        """
        Parse Primer3 output file and extract comprehensive primer information
        解析Primer3输出文件并提取全面的引物信息
        
        Args:
            input_file (str): Path to Primer3 output file / Primer3输出文件路径
            output_file (str): Path to output CSV file (optional) / 输出CSV文件路径（可选）
            
        Returns:
            pd.DataFrame: Parsed primer data with dimer parameters / 包含二聚体参数的解析引物数据
        """
        self.primers = []
        
        print(f"Parsing Primer3 output file: {input_file}")
        
        with open(input_file, 'r') as f:
            content = f.read()
        
        # Split into individual sequence results / 分割为单个序列结果
        sequences = content.split('=\n')
        
        for seq_block in sequences:
            if not seq_block.strip():
                continue
                
            self._parse_sequence_block(seq_block)
        
        # Convert to DataFrame / 转换为DataFrame
        df = pd.DataFrame(self.primers)
        
        print(f"Successfully parsed {len(df)} primer pairs")
        
        if output_file:
            df.to_csv(output_file, index=False)
            print(f"Results saved to {output_file}")
            
        return df
    
    def _parse_sequence_block(self, block):
        """
        Parse a single sequence block from Primer3 output / 解析Primer3输出中的单个序列块
        
        Args:
            block (str): Text block for one sequence / 单个序列的文本块
        """
        lines = block.strip().split('\n')
        sequence_data = {}
        
        # Parse sequence information / 解析序列信息
        for line in lines:
            if '=' in line:
                key, value = line.split('=', 1)
                sequence_data[key] = value
        
        # Extract sequence ID / 提取序列ID
        sequence_id = sequence_data.get('SEQUENCE_ID', 'unknown')
        
        # Find number of primer pairs returned / 查找返回的引物对数量
        num_pairs = int(sequence_data.get('PRIMER_PAIR_NUM_RETURNED', 0))
        
        if num_pairs == 0:
            print(f"Warning: No primer pairs found for sequence: {sequence_id}")
            return
            
        # Extract primer pairs / 提取引物对
        for i in range(num_pairs):
            primer_pair = self._extract_primer_pair(sequence_data, i, sequence_id)
            if primer_pair:
                self.primers.append(primer_pair)
    
    def _extract_primer_pair(self, data, pair_index, sequence_id):
        """
        Extract comprehensive information for a specific primer pair including dimer analysis
        提取特定引物对的全面信息，包括二聚体分析
        
        Args:
            data (dict): Parsed sequence data / 解析的序列数据
            pair_index (int): Index of primer pair (0, 1, 2...) / 引物对索引
            sequence_id (str): Sequence identifier / 序列标识符
            
        Returns:
            dict: Comprehensive primer pair information / 全面的引物对信息
        """
        try:
            # Extract basic left primer information / 提取基本左引物信息
            left_seq = data.get(f'PRIMER_LEFT_{pair_index}_SEQUENCE', '')
            left_tm = self._safe_float(data.get(f'PRIMER_LEFT_{pair_index}_TM', 0))
            left_gc = self._safe_float(data.get(f'PRIMER_LEFT_{pair_index}_GC_PERCENT', 0))
            left_pos = data.get(f'PRIMER_LEFT_{pair_index}', '0,0').split(',')
            left_start = int(left_pos[0]) if len(left_pos) > 0 else 0
            left_len = int(left_pos[1]) if len(left_pos) > 1 else len(left_seq)
            
            # Extract advanced left primer parameters / 提取高级左引物参数
            left_self_any = self._safe_float(data.get(f'PRIMER_LEFT_{pair_index}_SELF_ANY', 0))
            left_self_end = self._safe_float(data.get(f'PRIMER_LEFT_{pair_index}_SELF_END', 0))
            left_hairpin = self._safe_float(data.get(f'PRIMER_LEFT_{pair_index}_HAIRPIN', 0))
            left_end_stability = self._safe_float(data.get(f'PRIMER_LEFT_{pair_index}_END_STABILITY', 0))
            
            # Extract basic right primer information / 提取基本右引物信息
            right_seq = data.get(f'PRIMER_RIGHT_{pair_index}_SEQUENCE', '')
            right_tm = self._safe_float(data.get(f'PRIMER_RIGHT_{pair_index}_TM', 0))
            right_gc = self._safe_float(data.get(f'PRIMER_RIGHT_{pair_index}_GC_PERCENT', 0))
            right_pos = data.get(f'PRIMER_RIGHT_{pair_index}', '0,0').split(',')
            right_start = int(right_pos[0]) if len(right_pos) > 0 else 0
            right_len = int(right_pos[1]) if len(right_pos) > 1 else len(right_seq)
            
            # Extract advanced right primer parameters / 提取高级右引物参数
            right_self_any = self._safe_float(data.get(f'PRIMER_RIGHT_{pair_index}_SELF_ANY', 0))
            right_self_end = self._safe_float(data.get(f'PRIMER_RIGHT_{pair_index}_SELF_END', 0))
            right_hairpin = self._safe_float(data.get(f'PRIMER_RIGHT_{pair_index}_HAIRPIN', 0))
            right_end_stability = self._safe_float(data.get(f'PRIMER_RIGHT_{pair_index}_END_STABILITY', 0))
            
            # Extract primer pair complementarity / 提取引物对互补性
            pair_compl_any = self._safe_float(data.get(f'PRIMER_PAIR_{pair_index}_COMPL_ANY', 0))
            pair_compl_end = self._safe_float(data.get(f'PRIMER_PAIR_{pair_index}_COMPL_END', 0))
            
            # Extract product information / 提取产物信息
            product_size = int(data.get(f'PRIMER_PAIR_{pair_index}_PRODUCT_SIZE', 0))
            product_tm = self._safe_float(data.get(f'PRIMER_PAIR_{pair_index}_PRODUCT_TM', 0))
            
            # Calculate additional derived parameters / 计算额外的衍生参数
            tm_diff = abs(left_tm - right_tm)
            gc_diff = abs(left_gc - right_gc)
            len_diff = abs(left_len - right_len)
            
            primer_pair = {
                # Basic identification / 基本标识
                'ID': f"{sequence_id}_pair_{pair_index}",
                'Gene': sequence_id,
                'Pair_Index': pair_index,
                
                # Left primer basic info / 左引物基本信息
                'Left_Primer': left_seq,
                'Left_Tm': left_tm,
                'Left_GC%': left_gc,
                'Left_Start': left_start,
                'Left_Length': left_len,
                
                # Left primer dimer analysis / 左引物二聚体分析
                'Left_Self_Any': left_self_any,
                'Left_Self_End': left_self_end,
                'Left_Hairpin': left_hairpin,
                'Left_End_Stability': left_end_stability,
                
                # Right primer basic info / 右引物基本信息
                'Right_Primer': right_seq,
                'Right_Tm': right_tm,
                'Right_GC%': right_gc,
                'Right_Start': right_start,
                'Right_Length': right_len,
                
                # Right primer dimer analysis / 右引物二聚体分析
                'Right_Self_Any': right_self_any,
                'Right_Self_End': right_self_end,
                'Right_Hairpin': right_hairpin,
                'Right_End_Stability': right_end_stability,
                
                # Primer pair analysis / 引物对分析
                'Pair_Compl_Any': pair_compl_any,
                'Pair_Compl_End': pair_compl_end,
                'Product_Size': product_size,
                'Product_Tm': product_tm,
                
                # Derived parameters / 衍生参数
                'Tm_Diff': tm_diff,
                'GC_Diff': gc_diff,
                'Length_Diff': len_diff
            }
            
            return primer_pair
            
        except (ValueError, IndexError) as e:
            print(f"Error parsing primer pair {pair_index} for {sequence_id}: {e}")
            return None
    
    def _safe_float(self, value, default=0.0):
        """
        Safely convert value to float / 安全地将值转换为浮点数
        
        Args:
            value: Value to convert / 要转换的值
            default: Default value if conversion fails / 转换失败时的默认值
            
        Returns:
            float: Converted value or default / 转换后的值或默认值
        """
        try:
            return float(value) if value else default
        except (ValueError, TypeError):
            return default


def main():
    """
    Enhanced command line interface / 增强版命令行接口
    """
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Enhanced Primer3 output parser with comprehensive dimer analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例 / Examples:
  python primer3_parser.py primer3_output.txt parsed_primers.csv
  
功能特点 / Features:
  - 完整的引物参数提取 / Complete primer parameter extraction
  - 全面的二聚体分析 / Comprehensive dimer analysis
  - 质量指标计算 / Quality metrics calculation
  - 强健的错误处理 / Robust error handling
        """
    )
    
    parser.add_argument("input_file", help="Input Primer3 output file / 输入Primer3输出文件")
    parser.add_argument("output_file", help="Output CSV file / 输出CSV文件")
    parser.add_argument("--verbose", "-v", action="store_true", 
                       help="Enable verbose logging / 启用详细日志")
    
    args = parser.parse_args()
    
    # Auto replace spaces with underscores in file paths / 自动将路径中的空格替换为下划线
    input_file = args.input_file.replace(" ", "_")
    output_file = args.output_file.replace(" ", "_")
    
    # Parse the file / 解析文件
    primer_parser = Primer3Parser()
    df = primer_parser.parse_file(input_file, output_file)
    
    # Print summary / 打印摘要
    print(f"\n✅ Parsing completed successfully!")
    print(f"📊 Parsed {len(df)} primer pairs from {input_file}")
    print(f"📁 Results saved to {output_file}")
    
    if len(df) > 0:
        print(f"\n📋 Column summary:")
        basic_cols = [c for c in df.columns if not any(x in c for x in ['Self', 'Hairpin', 'Compl'])]
        dimer_cols = [c for c in df.columns if any(x in c for x in ['Self', 'Hairpin', 'Compl'])]
        print(f"   - Basic parameters: {len(basic_cols)} columns")
        print(f"   - Dimer parameters: {len(dimer_cols)} columns")
        
        # Show first few results / 显示前几个结果
        print(f"\n🔍 Sample results (first 3 rows):")
        sample_cols = ['ID', 'Left_Primer', 'Right_Primer', 'Product_Size', 'Left_Self_Any', 'Right_Self_Any']
        available_cols = [c for c in sample_cols if c in df.columns]
        if available_cols:
            print(df[available_cols].head(3).to_string(index=False))


if __name__ == "__main__":
    main() 