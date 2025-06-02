#!/usr/bin/env python3
"""
测试二聚体参数解析修复 / Test Dimer Parameter Parsing Fix

验证Primer3解析器和质量评分器的兼容性 / Verify compatibility between Primer3 parser and quality ranker
"""

import sys
import tempfile
from pathlib import Path

# Add toolkit to path
toolkit_path = Path(__file__).parent / "primer_design_toolkit"
sys.path.insert(0, str(toolkit_path))

from primer3_parser import Primer3Parser
from quality_ranker import PrimerQualityRanker


def create_mock_primer3_output():
    """创建模拟的Primer3输出 / Create mock Primer3 output"""
    
    mock_output = """SEQUENCE_ID=test_gene_1
 