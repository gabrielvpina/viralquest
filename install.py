# install.py
# This script installs the ViralQuest pipeline and its dependencies.
import os
import sys
import subprocess
import urllib.request
import tarfile
import shutil
from pathlib import Path


def run_command(cmd, check=True, shell=True):
    """Run a shell command and handle errors"""
    print(f"Running: {cmd}")
    try:
        result = subprocess.run(cmd, shell=shell, check=check, 
                              capture_output=True, text=True)
        if result.stdout:
            print(result.stdout)
        return result
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {cmd}")
        print(f"Error output: {e.stderr}")
        if check:
            sys.exit(1)
        return e


def check_conda():
    """Check if conda is available"""
    try:
        result = subprocess.run("conda --version", shell=True, 
                              capture_output=True, text=True)
        return result.returncode == 0
    except:
        return False


def get_conda_prefix():
    """Get the conda environment prefix"""
    try:
        result = subprocess.run("echo $CONDA_PREFIX", shell=True, 
                              capture_output=True, text=True)
        conda_prefix = result.stdout.strip()
        if conda_prefix and conda_prefix != "$CONDA_PREFIX":
            return conda_prefix
        
        # Alternative method
        result = subprocess.run("conda info --base", shell=True, 
                              capture_output=True, text=True)
        return result.stdout.strip()
    except:
        return None


def install_conda_packages():
    """Install required conda packages"""
    print("\n=== Installing conda packages ===")
    
    if not check_conda():
        print("Error: conda is not available. Please install conda/miniconda first.")
        sys.exit(1)
    
    # Install bioconda packages
    cmd = "conda install -c bioconda cap3 blast -y"
    run_command(cmd)


def install_diamond():
    """Download and install Diamond aligner"""
    print("\n=== Installing Diamond aligner ===")
    
    conda_prefix = get_conda_prefix()
    if not conda_prefix:
        print("Error: Could not determine conda prefix")
        sys.exit(1)
    
    bin_dir = os.path.join(conda_prefix, "bin")
    if not os.path.exists(bin_dir):
        print(f"Error: Conda bin directory not found: {bin_dir}")
        sys.exit(1)
    
    # Download Diamond
    diamond_url = "https://github.com/bbuchfink/diamond/releases/download/v2.1.12/diamond-linux64.tar.gz"
    diamond_tar = "diamond-linux64.tar.gz"
    
    print(f"Downloading Diamond from {diamond_url}")
    try:
        urllib.request.urlretrieve(diamond_url, diamond_tar)
    except Exception as e:
        print(f"Error downloading Diamond: {e}")
        sys.exit(1)
    
    # Extract Diamond
    print("Extracting Diamond...")
    try:
        with tarfile.open(diamond_tar, "r:gz") as tar:
            tar.extractall()
    except Exception as e:
        print(f"Error extracting Diamond: {e}")
        sys.exit(1)
    
    # Make executable and copy to bin
    diamond_binary = "./diamond"
    if os.path.exists(diamond_binary):
        os.chmod(diamond_binary, 0o755)
        destination = os.path.join(bin_dir, "diamond")
        shutil.copy2(diamond_binary, destination)
        print(f"Diamond installed to {destination}")
        
        # Cleanup
        os.remove(diamond_binary)
        os.remove(diamond_tar)
    else:
        print("Error: Diamond binary not found after extraction")
        sys.exit(1)


def install_pip_requirements():
    """Install pip packages from requirements.txt"""
    print("\n=== Installing pip packages ===")
    
    if not os.path.exists("requirements.txt"):
        print("Warning: requirements.txt not found. Skipping pip installation.")
        return
    
    cmd = "pip install -r requirements.txt"
    run_command(cmd)


def setup_viralquest_executable():
    """Make viralquest.py executable and copy to conda bin"""
    print("\n=== Setting up ViralQuest executable ===")
    
    if not os.path.exists("viralquest.py"):
        print("Error: viralquest.py not found in current directory")
        sys.exit(1)
    
    conda_prefix = get_conda_prefix()
    if not conda_prefix:
        print("Error: Could not determine conda prefix")
        sys.exit(1)
    
    bin_dir = os.path.join(conda_prefix, "bin")
    
    # Make executable
    os.chmod("viralquest.py", 0o755)
    print("Made viralquest.py executable")
    
    # Copy to conda bin
    destination = os.path.join(bin_dir, "viralquest.py")
    shutil.copy2("viralquest.py", destination)
    print(f"Copied viralquest.py to {destination}")


def test_installation():
    """Test if the installation was successful"""
    print("\n=== Testing installation ===")
    
    # test viralquest.py
    print("Testing viralquest.py...")
    result = run_command("python viralquest.py --help", check=False)
    if result.returncode == 0:
        print("✓ viralquest.py works with python")
    else:
        print("✗ viralquest.py failed with python")
    
    # test executable
    print("Testing viralquest.py as executable...")
    result = run_command("viralquest.py --help", check=False)
    if result.returncode == 0:
        print("✓ viralquest.py works as executable")
    else:
        print("✗ viralquest.py failed as executable")
    
    # test dependencies
    print("Testing dependencies...")
    dependencies = ["cap3", "blastn", "diamond"]
    for dep in dependencies:
        result = run_command(f"which {dep}", check=False)
        if result.returncode == 0:
            print(f"✓ {dep} found")
        else:
            print(f"✗ {dep} not found")


def main():
    """Main setup function"""
    print("ViralQuest Setup Script")
    print("=" * 50)
    
    if not os.path.exists("viralquest.py"):
        print("Error: This script should be run from the directory containing viralquest.py")
        sys.exit(1)
    
    try:
        install_conda_packages()
        
        install_diamond()
        
        install_pip_requirements()
        
        setup_viralquest_executable()
        
        test_installation()
        
        print("\n" + "=" * 50)
        print("Setup completed successfully!")
        print("You can now run: viralquest.py --help")
        
    except KeyboardInterrupt:
        print("\nSetup interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\nUnexpected error during setup: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()