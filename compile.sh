find . -name \*.pyc -exec rm {} \;
python2 -c "import compileall; compileall.compile_dir('.', force=True)"

