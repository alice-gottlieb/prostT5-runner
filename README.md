# Installing ProstT5 and Foldseek on Hoffman2 #
Important: When installing packages, set `uv pip install --only-binary=:all: [package name]`. This is necessary for Hoffman to install and run the packages correctly

1.  Install uv
2. Activate uv venv provided
3. Install foldseek 
`wget https://mmseqs.com/foldseek/foldseek-linux-gpu.tar.gz; tar xvfz foldseek-linux-gpu.tar.gz; export PATH=$(pwd)/foldseek/bin/:$PATH`
