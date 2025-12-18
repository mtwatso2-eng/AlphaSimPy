#!/bin/bash
# Simple installation script for AlphaSimPy

echo "AlphaSimPy Installation"
echo "======================"

# Check if we're in the right directory
if [ ! -f "AlphaSimPy.py" ]; then
    echo "Error: AlphaSimPy.py not found. Please run this script from the AlphaSimPy directory."
    exit 1
fi

# Check for Python
if ! command -v python3 &> /dev/null; then
    echo "Error: Python 3 is required but not installed."
    exit 1
fi

# Check for pip
if ! command -v pip3 &> /dev/null; then
    echo "Error: pip3 is required but not installed."
    exit 1
fi

# Install Python dependencies
echo "Installing Python dependencies..."
pip3 install -r requirements.txt

# Test the installation
echo "Testing the installation..."
python3 -c "from AlphaSimPy import runMacs; print('✓ AlphaSimPy imported successfully')"

echo ""
echo "Installation completed successfully!"
echo "You can now use AlphaSimPy by running:"
echo "  python3 test_alphasimpy.py"
echo "  or"
echo "  jupyter notebook AlphaSimPy.ipynb"
