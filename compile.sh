find . -name \*.pyc -exec rm {} \;
python -c "import compileall; compileall.compile_dir('.', force=True)"

