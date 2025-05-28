#!/usr/bin/env python3
"""
Primer3 Result Parser
Parse Primer3 output files and extract primer information
"""

import pandas as pd
import re
from pathlib import Path


class Primer3Parser:
    """
    Parser for Primer3 output files
    Extracts primer pair information and converts to structured format
    """
    
    def __init__(self):
        self.primers = []
        
    def parse_file(self, input_file, output_file=None):
        """
        Parse Primer3 output file and extract primer information
        
        Args:
            input_file (str): Path to Primer3 output file
            output_file (str): Path to output CSV file (optional)
            
        Returns:
            pd.DataFrame: Parsed primer data
        """
        self.primers = []
        
        with open(input_file, 'r') as f:
            content = f.read()
        
        # Split into individual sequence results
        sequences = content.split('=\n')
        
        for seq_block in sequences:
            if not seq_block.strip():
                continue
                
            self._parse_sequence_block(seq_block)
        
        # Convert to DataFrame
        df = pd.DataFrame(self.primers)
        
        if output_file:
            df.to_csv(output_file, index=False)
            
        return df
    
    def _parse_sequence_block(self, block):
        """
        Parse a single sequence block from Primer3 output
        
        Args:
            block (str): Text block for one sequence
        """
        lines = block.strip().split('\n')
        sequence_data = {}
        
        # Parse sequence information
        for line in lines:
            if '=' in line:
                key, value = line.split('=', 1)
                sequence_data[key] = value
        
        # Extract sequence ID
        sequence_id = sequence_data.get('SEQUENCE_ID', 'unknown')
        
        # Find number of primer pairs returned
        num_pairs = int(sequence_data.get('PRIMER_PAIR_NUM_RETURNED', 0))
        
        if num_pairs == 0:
            return
            
        # Extract primer pairs
        for i in range(num_pairs):
            primer_pair = self._extract_primer_pair(sequence_data, i, sequence_id)
            if primer_pair:
                self.primers.append(primer_pair)
    
    def _extract_primer_pair(self, data, pair_index, sequence_id):
        """
        Extract information for a specific primer pair
        
        Args:
            data (dict): Parsed sequence data
            pair_index (int): Index of primer pair (0, 1, 2...)
            sequence_id (str): Sequence identifier
            
        Returns:
            dict: Primer pair information
        """
        try:
            # Extract left primer
            left_seq = data.get(f'PRIMER_LEFT_{pair_index}_SEQUENCE', '')
            left_tm = float(data.get(f'PRIMER_LEFT_{pair_index}_TM', 0))
            left_gc = float(data.get(f'PRIMER_LEFT_{pair_index}_GC_PERCENT', 0))
            left_pos = data.get(f'PRIMER_LEFT_{pair_index}', '0,0').split(',')
            left_start = int(left_pos[0]) if len(left_pos) > 0 else 0
            left_len = int(left_pos[1]) if len(left_pos) > 1 else len(left_seq)
            
            # Extract right primer
            right_seq = data.get(f'PRIMER_RIGHT_{pair_index}_SEQUENCE', '')
            right_tm = float(data.get(f'PRIMER_RIGHT_{pair_index}_TM', 0))
            right_gc = float(data.get(f'PRIMER_RIGHT_{pair_index}_GC_PERCENT', 0))
            right_pos = data.get(f'PRIMER_RIGHT_{pair_index}', '0,0').split(',')
            right_start = int(right_pos[0]) if len(right_pos) > 0 else 0
            right_len = int(right_pos[1]) if len(right_pos) > 1 else len(right_seq)
            
            # Extract product information
            product_size = int(data.get(f'PRIMER_PAIR_{pair_index}_PRODUCT_SIZE', 0))
            
            primer_pair = {
                'ID': f"{sequence_id}_pair_{pair_index}",
                'Gene': sequence_id,
                'Pair_Index': pair_index,
                'Left_Primer': left_seq,
                'Left_Tm': left_tm,
                'Left_GC': left_gc,
                'Left_Start': left_start,
                'Left_Length': left_len,
                'Right_Primer': right_seq,
                'Right_Tm': right_tm,
                'Right_GC': right_gc,
                'Right_Start': right_start,
                'Right_Length': right_len,
                'Product_Size': product_size
            }
            
            return primer_pair
            
        except (ValueError, IndexError) as e:
            print(f"Error parsing primer pair {pair_index} for {sequence_id}: {e}")
            return None


def main():
    """
    Command line interface for Primer3Parser
    """
    import argparse
    
    parser = argparse.ArgumentParser(description="Parse Primer3 output files")
    parser.add_argument("input_file", help="Input Primer3 output file")
    parser.add_argument("output_file", help="Output CSV file")
    
    args = parser.parse_args()
    
    # Parse the file
    p3_parser = Primer3Parser()
    df = p3_parser.parse_file(args.input_file, args.output_file)
    
    print(f"Parsed {len(df)} primer pairs from {args.input_file}")
    print(f"Results saved to {args.output_file}")


if __name__ == "__main__":
    main() 