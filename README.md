# Utilities for metabolic modeling 

`fluxpy` is a toolkit for simple but laborious tasks in metabolic modeling analysis. 


## Install 

<!-- in the .toml file we keep the dependencies and that's it.  -->


If you wish just to use `fluxpy`:

```
pip install fluxpy
```

----

If you wish to contribute/develop, git clone this repo, fire a new branch, and after you add your changes, you may run: 

```
python -m build

pip install .
```

You can then open a PR to ask your code to be merged on the library! 


## ReadTheDocs

To test whether the documentation builder will perform online as we would like to, we may run locally the following command:

```bash
cd docs/
sphinx-build -b html -d _build/doctrees -D language=en . _build/html -v
```


