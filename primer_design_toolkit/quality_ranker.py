#!/usr/bin/env python3
"""
Primer Quality Ranker - Enhanced Version
引物质量评分器 - 增强版 / Primer Quality Ranker - Enhanced Version

Comprehensive primer quality evaluation with parallel processing
基于经过验证的质量指标和评分算法的全面引物质量评估，支持并行处理
"""

import pandas as pd
import numpy as np
from pathlib import Path
import concurrent.futures
from typing import Dict, List, Tuple
import logging


class PrimerQualityRanker:
    """
    引物质量评估和排序系统 - 增强版 / Enhanced primer quality evaluation and ranking system
    
    使用全面的质量指标进行并行处理 / Uses comprehensive quality metrics with parallel processing
    基于经过验证的10维评分体系 / Based on proven 10-dimensional scoring system
    """
    
    def __init__(self, threads=4):
        """
        初始化质量评分器，包含全面的参数配置 / Initialize the quality ranker with comprehensive parameters
        
        Args:
            threads (int): 并行处理线程数 / Number of threads for parallel processing
        """
        self.threads = threads
        
        # 全面的质量评估参数（基于经过验证的指标）/ Comprehensive quality evaluation parameters (based on proven metrics)
        self.IDEAL_PARAMS = {
            'product_size': {'ideal': 100, 'weight': 0.10, 'range': (80, 200)},    # 产物大小 / Product size
            'tm': {'ideal': 60, 'weight': 0.15, 'range': (55, 65)},                # 熔解温度 / Melting temperature
            'gc': {'ideal': 50, 'weight': 0.15, 'range': (40, 60)},                # GC含量 / GC content
            'tm_diff': {'ideal': 0, 'weight': 0.10, 'max_diff': 5},                # Tm差异 / Tm difference
            'len': {'ideal': 20, 'weight': 0.10, 'range': (18, 25)},               # 引物长度 / Primer length
            'gc_diff': {'ideal': 0, 'weight': 0.05, 'max_diff': 10},               # GC差异 / GC difference
            # 关键的二聚体和结构参数 / Critical dimer and structure parameters
            'self_any': {'ideal': 0, 'weight': 0.10, 'max': 8.0},                  # 自身二聚体 / Self dimer
            'self_end': {'ideal': 0, 'weight': 0.10, 'max': 3.0},                  # 末端二聚体 / End dimer
            'hairpin': {'ideal': 0, 'weight': 0.05, 'max': 40.0},                  # 发夹结构 / Hairpin structure
            'pair_compl': {'ideal': 0, 'weight': 0.10, 'max': 8.0}                 # 引物对互补性 / Primer pair complementarity
        }
        
        # 设置日志记录 / Setup logging
        self.logger = logging.getLogger(__name__)
        
    def rank_primers(self, input_file, output_file=None):
        """
        基于全面质量指标和并行处理对引物对进行排序 / Rank primer pairs based on comprehensive quality metrics with parallel processing
        
        Args:
            input_file (str): 解析后的引物CSV文件路径 / Path to parsed primer CSV file
            output_file (str): 输出排序后CSV文件路径（可选）/ Path to output ranked CSV file (optional)
            
        Returns:
            pd.DataFrame: 包含质量得分的排序引物数据 / Ranked primer data with quality scores
        """
        self.logger.info(f"Loading primer data from {input_file}")
        
        # 加载引物数据 / Load primer data
        df = pd.read_csv(input_file)
        
        if df.empty:
            self.logger.warning("No primer data found in input file")
            return df
        
        self.logger.info(f"Evaluating quality for {len(df)} primer pairs using {self.threads} threads")
        
        # 使用并行处理计算质量得分 / Calculate quality scores with parallel processing
        df = self._calculate_quality_scores_parallel(df)
        
        # 添加分组排序和全局排序 / Add group ranking and global ranking
        df = self._add_rankings(df)
        
        # 添加质量等级 / Add quality grades
        df = self._add_quality_grades(df)
        
        # 重新排列列顺序以提高可读性 / Reorder columns for better readability
        df = self._reorder_columns(df)
        
        # 按质量得分排序（降序）/ Sort by quality score (descending)
        df = df.sort_values('Quality_Score', ascending=False).reset_index(drop=True)
        
        if output_file:
            df.to_csv(output_file, index=False)
            self.logger.info(f"Results saved to {output_file}")
            
        # 打印统计摘要 / Print summary statistics
        self._print_quality_summary(df)
        
        return df
    
    def _calculate_quality_scores_parallel(self, df):
        """
        使用并行处理计算质量得分 / Calculate quality scores using parallel processing
        
        Args:
            df (pd.DataFrame): 引物数据 / Primer data
            
        Returns:
            pd.DataFrame: 添加了质量指标的数据 / Data with added quality metrics
        """
        # 将数据框分割成块进行并行处理 / Split dataframe into chunks for parallel processing
        chunk_size = max(1, len(df) // self.threads)
        chunks = [df.iloc[i:i + chunk_size] for i in range(0, len(df), chunk_size)]
        
        # 并行处理块 / Process chunks in parallel
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.threads) as executor:
            future_to_chunk = {
                executor.submit(self._calculate_quality_scores_chunk, chunk): i 
                for i, chunk in enumerate(chunks)
            }
            
            processed_chunks = []
            for future in concurrent.futures.as_completed(future_to_chunk):
                chunk_idx = future_to_chunk[future]
                try:
                    processed_chunk = future.result()
                    processed_chunks.append((chunk_idx, processed_chunk))
                except Exception as e:
                    self.logger.error(f"Error processing chunk {chunk_idx}: {e}")
                    raise
        
        # 按原始顺序合并结果 / Combine results in original order
        processed_chunks.sort(key=lambda x: x[0])
        result_df = pd.concat([chunk for _, chunk in processed_chunks], ignore_index=True)
        
        return result_df
    
    def _calculate_quality_scores_chunk(self, df_chunk):
        """
        计算数据块的质量得分 / Calculate quality scores for a chunk of data
        
        Args:
            df_chunk (pd.DataFrame): 引物数据块 / Chunk of primer data
            
        Returns:
            pd.DataFrame: 添加了质量指标的数据块 / Chunk with added quality metrics
        """
        # 对每行应用质量评分 / Apply quality scoring to each row
        df_chunk = df_chunk.copy()
        df_chunk['Quality_Score'] = df_chunk.apply(self._calculate_single_quality_score, axis=1)
        
        # 如果可用，添加二聚体特定得分 / Add dimer-specific score if available
        if self._has_dimer_columns(df_chunk):
            df_chunk['Dimer_Score'] = df_chunk.apply(self._calculate_dimer_score, axis=1)
        
        return df_chunk
    
    def _calculate_single_quality_score(self, row):
        """
        计算单个引物对的综合质量得分 / Calculate comprehensive quality score for a single primer pair
        基于经过验证的10维评估算法 / Based on proven 10-dimensional evaluation algorithm
        
        Args:
            row (pd.Series): 单个引物对数据 / Single primer pair data
            
        Returns:
            float: 质量得分 (0-100) / Quality score (0-100)
        """
        scores = []
        
        # 1. 产物大小得分 / Product size score
        product_size = row['Product_Size']
        ideal_size = self.IDEAL_PARAMS['product_size']['ideal']
        size_range = self.IDEAL_PARAMS['product_size']['range']
        
        if size_range[0] <= product_size <= size_range[1]:
            size_score = 1.0 - min(abs(product_size - ideal_size), 50) / 50
        else:
            size_score = 0.3
        scores.append(size_score * self.IDEAL_PARAMS['product_size']['weight'])
        
        # 2. Tm得分（左右引物平均）/ Tm score (average of left and right)
        left_tm = row['Left_Tm']
        right_tm = row['Right_Tm']
        ideal_tm = self.IDEAL_PARAMS['tm']['ideal']
        tm_range = self.IDEAL_PARAMS['tm']['range']
        
        left_tm_score = 1.0 - abs(left_tm - ideal_tm) / 10 if tm_range[0] <= left_tm <= tm_range[1] else 0.3
        right_tm_score = 1.0 - abs(right_tm - ideal_tm) / 10 if tm_range[0] <= right_tm <= tm_range[1] else 0.3
        tm_score = (left_tm_score + right_tm_score) / 2
        scores.append(tm_score * self.IDEAL_PARAMS['tm']['weight'])
        
        # 3. Tm差异得分 / Tm difference score
        tm_diff = abs(left_tm - right_tm)
        max_tm_diff = self.IDEAL_PARAMS['tm_diff']['max_diff']
        tm_diff_score = 1.0 - min(tm_diff, max_tm_diff) / max_tm_diff
        scores.append(tm_diff_score * self.IDEAL_PARAMS['tm_diff']['weight'])
        
        # 4. GC含量得分 / GC content score
        left_gc = row.get('Left_GC%', row.get('Left_GC', 0))
        right_gc = row.get('Right_GC%', row.get('Right_GC', 0))
        ideal_gc = self.IDEAL_PARAMS['gc']['ideal']
        gc_range = self.IDEAL_PARAMS['gc']['range']
        
        left_gc_score = 1.0 - abs(left_gc - ideal_gc) / 20 if gc_range[0] <= left_gc <= gc_range[1] else 0.3
        right_gc_score = 1.0 - abs(right_gc - ideal_gc) / 20 if gc_range[0] <= right_gc <= gc_range[1] else 0.3
        gc_score = (left_gc_score + right_gc_score) / 2
        scores.append(gc_score * self.IDEAL_PARAMS['gc']['weight'])
        
        # 5. GC差异得分 / GC difference score
        gc_diff = abs(left_gc - right_gc)
        max_gc_diff = self.IDEAL_PARAMS['gc_diff']['max_diff']
        gc_diff_score = 1.0 - min(gc_diff, max_gc_diff) / max_gc_diff
        scores.append(gc_diff_score * self.IDEAL_PARAMS['gc_diff']['weight'])
        
        # 6. 引物长度得分 / Primer length score
        left_len = row.get('Left_Len', row.get('Left_Length', 0))
        right_len = row.get('Right_Len', row.get('Right_Length', 0))
        ideal_len = self.IDEAL_PARAMS['len']['ideal']
        len_range = self.IDEAL_PARAMS['len']['range']
        
        left_len_score = 1.0 - abs(left_len - ideal_len) / 10 if len_range[0] <= left_len <= len_range[1] else 0.3
        right_len_score = 1.0 - abs(right_len - ideal_len) / 10 if len_range[0] <= right_len <= len_range[1] else 0.3
        len_score = (left_len_score + right_len_score) / 2
        scores.append(len_score * self.IDEAL_PARAMS['len']['weight'])
        
        # 7. 自身二聚体得分 (Self_Any) / Self-dimer score (Self_Any)
        if 'Left_Self_Any' in row and 'Right_Self_Any' in row:
            left_self_any = row['Left_Self_Any']
            right_self_any = row['Right_Self_Any']
            max_self_any = self.IDEAL_PARAMS['self_any']['max']
            
            left_self_any_score = 1.0 - min(left_self_any, max_self_any) / max_self_any
            right_self_any_score = 1.0 - min(right_self_any, max_self_any) / max_self_any
            self_any_score = (left_self_any_score + right_self_any_score) / 2
            scores.append(self_any_score * self.IDEAL_PARAMS['self_any']['weight'])
        
        # 8. 末端二聚体得分 (Self_End) / End-dimer score (Self_End)
        if 'Left_Self_End' in row and 'Right_Self_End' in row:
            left_self_end = row['Left_Self_End']
            right_self_end = row['Right_Self_End']
            max_self_end = self.IDEAL_PARAMS['self_end']['max']
            
            left_self_end_score = 1.0 - min(left_self_end, max_self_end) / max_self_end
            right_self_end_score = 1.0 - min(right_self_end, max_self_end) / max_self_end
            self_end_score = (left_self_end_score + right_self_end_score) / 2
            scores.append(self_end_score * self.IDEAL_PARAMS['self_end']['weight'])
        
        # 9. 发夹结构得分 / Hairpin structure score
        if 'Left_Hairpin' in row and 'Right_Hairpin' in row:
            left_hairpin = row['Left_Hairpin']
            right_hairpin = row['Right_Hairpin']
            max_hairpin = self.IDEAL_PARAMS['hairpin']['max']
            
            left_hairpin_score = 1.0 - min(left_hairpin, max_hairpin) / max_hairpin
            right_hairpin_score = 1.0 - min(right_hairpin, max_hairpin) / max_hairpin
            hairpin_score = (left_hairpin_score + right_hairpin_score) / 2
            scores.append(hairpin_score * self.IDEAL_PARAMS['hairpin']['weight'])
        
        # 10. 引物对互补性得分 / Primer pair complementarity score
        if 'Pair_Compl_Any' in row and 'Pair_Compl_End' in row:
            pair_compl_any = row['Pair_Compl_Any']
            pair_compl_end = row['Pair_Compl_End']
            max_pair_compl = self.IDEAL_PARAMS['pair_compl']['max']
            
            pair_compl_any_score = 1.0 - min(pair_compl_any, max_pair_compl) / max_pair_compl
            pair_compl_end_score = 1.0 - min(pair_compl_end * 2, max_pair_compl) / max_pair_compl
            pair_compl_score = (pair_compl_any_score * 0.4 + pair_compl_end_score * 0.6)
            scores.append(pair_compl_score * self.IDEAL_PARAMS['pair_compl']['weight'])
        
        # 总分 (0-100) / Total score (0-100)
        total_score = sum(scores) * 100
        return min(100, max(0, total_score))  # 确保得分在0-100范围内 / Ensure score is within 0-100 range
    
    def _calculate_dimer_score(self, row):
        """
        计算用于详细分析的特定二聚体得分 / Calculate specific dimer score for detailed analysis
        
        Args:
            row (pd.Series): 单个引物对数据 / Single primer pair data
            
        Returns:
            float: 二聚体得分 (0-100) / Dimer score (0-100)
        """
        if not self._has_dimer_columns(pd.DataFrame([row])):
            return 0
        
        left_self_any = row.get('Left_Self_Any', 0)
        right_self_any = row.get('Right_Self_Any', 0)
        pair_compl_any = row.get('Pair_Compl_Any', 0)
        
        # 基于经过验证的公式计算二聚体得分 / Calculate dimer score based on proven formula
        dimer_score = 100 * (1 - (left_self_any + right_self_any + pair_compl_any * 2) / 24)
        return min(100, max(0, dimer_score))
    
    def _has_dimer_columns(self, df):
        """检查数据框是否包含二聚体相关列 / Check if dataframe has dimer-related columns"""
        dimer_cols = ['Left_Self_Any', 'Right_Self_Any', 'Pair_Compl_Any']
        return all(col in df.columns for col in dimer_cols)
    
    def _add_rankings(self, df):
        """添加分组和全局排序 / Add group and global rankings"""
        # 按ID分组排序 / Group ranking by ID
        if 'ID' in df.columns:
            df['Rank_In_Group'] = df.groupby('ID')['Quality_Score'].rank(
                ascending=False, method='first'
            ).astype(int)
        
        # 全局排序 / Global ranking
        df['Global_Rank'] = df['Quality_Score'].rank(ascending=False, method='first').astype(int)
        
        return df
    
    def _add_quality_grades(self, df):
        """基于得分范围添加质量等级 / Add quality grades based on score ranges"""
        def get_grade(score):
            if score >= 90:
                return 'A+'
            elif score >= 85:
                return 'A'
            elif score >= 80:
                return 'B+'
            elif score >= 75:
                return 'B'
            elif score >= 70:
                return 'C+'
            elif score >= 65:
                return 'C'
            else:
                return 'D'
        
        df['Grade'] = df['Quality_Score'].apply(get_grade)
        return df
    
    def _reorder_columns(self, df):
        """重新排列列以提高可读性 / Reorder columns for better readability"""
        # 定义优先列顺序 / Define preferred column order
        priority_columns = [
            'ID', 'Primer_Pair_Index', 'Quality_Score', 'Grade', 
            'Global_Rank', 'Rank_In_Group', 'Product_Size',
            'Left_Seq', 'Left_Start', 'Left_Len', 'Left_Tm', 'Left_GC%',
            'Right_Seq', 'Right_Start', 'Right_Len', 'Right_Tm', 'Right_GC%'
        ]
        
        # 如果存在，添加二聚体列 / Add dimer columns if they exist
        dimer_columns = [
            'Left_Self_Any', 'Right_Self_Any', 'Left_Self_End', 'Right_Self_End',
            'Left_Hairpin', 'Right_Hairpin', 'Pair_Compl_Any', 'Pair_Compl_End',
            'Dimer_Score'
        ]
        
        # 构建最终列顺序 / Build final column order
        final_columns = []
        for col in priority_columns:
            if col in df.columns:
                final_columns.append(col)
        
        for col in dimer_columns:
            if col in df.columns:
                final_columns.append(col)
        
        # 添加任何剩余的列 / Add any remaining columns
        for col in df.columns:
            if col not in final_columns:
                final_columns.append(col)
        
        return df[final_columns]
    
    def _print_quality_summary(self, df):
        """打印综合质量摘要 / Print comprehensive quality summary"""
        total_primers = len(df)
        high_quality = len(df[df['Quality_Score'] >= 80])
        medium_quality = len(df[(df['Quality_Score'] >= 70) & (df['Quality_Score'] < 80)])
        low_quality = len(df[df['Quality_Score'] < 70])
        
        print(f"\nPrimer quality evaluation completed!")
        print(f"Total primers evaluated: {total_primers}")
        print(f"High quality primers (≥80): {high_quality} ({high_quality/total_primers*100:.1f}%)")
        print(f"Medium quality primers (70-79): {medium_quality} ({medium_quality/total_primers*100:.1f}%)")
        print(f"Low quality primers (<70): {low_quality} ({low_quality/total_primers*100:.1f}%)")
        
        if 'Dimer_Score' in df.columns:
            good_dimer = len(df[df['Dimer_Score'] >= 80])
            print(f"Good dimer performance (≥80): {good_dimer} ({good_dimer/total_primers*100:.1f}%)")
        
        # 显示前5名引物 / Show top 5 primers
        if len(df) > 0:
            print(f"\nTop 5 primers:")
            top_5 = df.head(5)[['ID', 'Quality_Score', 'Grade', 'Product_Size']].to_string(index=False)
            print(top_5)
    
    def get_quality_summary(self, df):
        """
        生成详细的质量摘要统计 / Generate detailed quality summary statistics
        
        Args:
            df (pd.DataFrame): 排序后的引物数据 / Ranked primer data
            
        Returns:
            dict: 综合质量摘要 / Comprehensive quality summary
        """
        if df.empty:
            return {}
        
        summary = {
            'total_primers': len(df),
            'high_quality': len(df[df['Quality_Score'] >= 80]),
            'medium_quality': len(df[(df['Quality_Score'] >= 60) & (df['Quality_Score'] < 80)]),
            'low_quality': len(df[df['Quality_Score'] < 60]),
            'avg_quality_score': df['Quality_Score'].mean(),
            'max_quality_score': df['Quality_Score'].max(),
            'min_quality_score': df['Quality_Score'].min(),
            'avg_tm_left': df['Left_Tm'].mean(),
            'avg_tm_right': df['Right_Tm'].mean(),
            'avg_product_size': df['Product_Size'].mean()
        }
        
        # 添加GC含量统计 / Add GC content stats
        if 'Left_GC%' in df.columns:
            summary['avg_gc_left'] = df['Left_GC%'].mean()
            summary['avg_gc_right'] = df['Right_GC%'].mean()
        
        # 如果可用，添加二聚体统计 / Add dimer stats if available
        if 'Dimer_Score' in df.columns:
            summary['avg_dimer_score'] = df['Dimer_Score'].mean()
            summary['good_dimer_count'] = len(df[df['Dimer_Score'] >= 80])
        
        return summary


def main():
    """
    命令行接口 / Command line interface
    """
    import argparse
    parser = argparse.ArgumentParser(description="Primer quality ranking")
    parser.add_argument("input_file", help="Input parsed primer CSV file")
    parser.add_argument("output_file", help="Output ranked CSV file")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads")
    args = parser.parse_args()

    # 自动将路径中的空格替换为下划线 / Auto replace spaces with underscores in file paths
    input_file = args.input_file.replace(" ", "_")
    output_file = args.output_file.replace(" ", "_")

    ranker = PrimerQualityRanker(threads=args.threads)
    ranker.rank_primers(input_file, output_file)


if __name__ == "__main__":
    main() 