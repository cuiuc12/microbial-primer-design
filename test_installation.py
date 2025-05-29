#!/usr/bin/env python3
"""
安装测试脚本 / Installation Test Script

测试微生物引物设计工具包的安装和依赖
Test the installation and dependencies of Microbial Primer Design Toolkit
"""

import sys
import subprocess
from pathlib import Path

def test_python_dependencies():
    """测试Python依赖 / Test Python dependencies"""
    print("🔍 Testing Python dependencies...")
    
    required_packages = ['pandas', 'biopython', 'numpy']
    missing_packages = []
    
    for package in required_packages:
        try:
            __import__(package)
            print(f"   ✅ {package}")
        except ImportError:
            print(f"   ❌ {package} - Missing")
            missing_packages.append(package)
    
    return len(missing_packages) == 0, missing_packages

def test_external_tools():
    """测试外部工具 / Test external tools"""
    print("\n🔍 Testing external tools...")
    
    tools = {
        'prokka': ['prokka', '--version'],
        'roary': ['roary', '--version'],
        'mafft': ['mafft', '--version'],
        'primer3': ['primer3_core', '--version']
    }
    
    missing_tools = []
    
    for tool_name, command in tools.items():
        try:
            result = subprocess.run(command, capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                print(f"   ✅ {tool_name}")
            else:
                print(f"   ❌ {tool_name} - Command failed")
                missing_tools.append(tool_name)
        except (subprocess.TimeoutExpired, FileNotFoundError):
            print(f"   ❌ {tool_name} - Not found")
            missing_tools.append(tool_name)
    
    return len(missing_tools) == 0, missing_tools

def test_toolkit_import():
    """测试工具包导入 / Test toolkit import"""
    print("\n🔍 Testing toolkit import...")
    
    try:
        from primer_design_toolkit import GenomePipeline, GenomeDownloader, Primer3Parser, PrimerQualityRanker
        print("   ✅ All toolkit modules imported successfully")
        return True
    except ImportError as e:
        print(f"   ❌ Import failed: {e}")
        return False

def test_file_structure():
    """测试文件结构 / Test file structure"""
    print("\n🔍 Testing file structure...")
    
    required_files = [
        'run_primer_design.py',
        'primer_design_toolkit/__init__.py',
        'primer_design_toolkit/primer_pipeline.py',
        'primer_design_toolkit/genome_downloader.py',
        'primer_design_toolkit/primer3_parser.py',
        'primer_design_toolkit/quality_ranker.py',
        'utils/__init__.py',
        'utils/download_genomes.py',
        'utils/parse_primer3.py',
        'utils/rank_primers.py',
        'requirements.txt',
        'environment.yml',
        'LICENSE',
        'README.md'
    ]
    
    missing_files = []
    
    for file_path in required_files:
        if Path(file_path).exists():
            print(f"   ✅ {file_path}")
        else:
            print(f"   ❌ {file_path} - Missing")
            missing_files.append(file_path)
    
    return len(missing_files) == 0, missing_files

def main():
    """主测试函数 / Main test function"""
    print("🧪 Microbial Primer Design Toolkit - Installation Test")
    print("=" * 60)
    
    all_tests_passed = True
    
    # 测试文件结构 / Test file structure
    files_ok, missing_files = test_file_structure()
    if not files_ok:
        all_tests_passed = False
    
    # 测试Python依赖 / Test Python dependencies
    deps_ok, missing_deps = test_python_dependencies()
    if not deps_ok:
        all_tests_passed = False
    
    # 测试工具包导入 / Test toolkit import
    import_ok = test_toolkit_import()
    if not import_ok:
        all_tests_passed = False
    
    # 测试外部工具 / Test external tools
    tools_ok, missing_tools = test_external_tools()
    if not tools_ok:
        all_tests_passed = False
    
    # 总结 / Summary
    print("\n" + "=" * 60)
    if all_tests_passed:
        print("🎉 All tests passed! Installation is complete and ready to use.")
    else:
        print("❌ Some tests failed. Please check the issues above.")
        
        if missing_files:
            print(f"\n📁 Missing files: {', '.join(missing_files)}")
        
        if missing_deps:
            print(f"\n🐍 Missing Python packages: {', '.join(missing_deps)}")
            print("   Install with: pip install " + " ".join(missing_deps))
        
        if missing_tools:
            print(f"\n🔧 Missing external tools: {', '.join(missing_tools)}")
            print("   Install with: conda install -c bioconda " + " ".join(missing_tools))
    
    print("=" * 60)
    return all_tests_passed

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1) 