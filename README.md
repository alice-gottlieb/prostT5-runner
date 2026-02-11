# Installing ProstT5 and Foldseek on Hoffman2 #
Important: When installing packages, set `uv pip install --only-binary=:all: [package name]`. This is necessary for Hoffman to install and run the packages correctly

Important 2: Run all python files with `uv run python [name of file].py`. This is the only way I've found so far to get uv to correctly apply dependencies.

1.  Install uv
2. Activate uv venv provided
3. Install foldseek 
`wget https://mmseqs.com/foldseek/foldseek-linux-gpu.tar.gz; tar xvfz foldseek-linux-gpu.tar.gz; export PATH=$(pwd)/foldseek/bin/:$PATH`
