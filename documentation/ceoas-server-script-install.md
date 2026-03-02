# Notes on getting the CEOAS server environment set up for the PIV pipline.

## Getting the scripts

1. Probably best to just clone my github repo here:

https://github.com/shindelr/scripts

2. Or, if that doesn't work, use these paths to access them from Robin's user

/home/server/pi/homes/shindelr/scripts/pivpipe
/home/server/pi/homes/shindelr/scripts/batchscript


## PyEnv stuff

To run these scripts, Python3.11 is recommended. CEOAS seems to have Python3.7 set up by default. The hacky way to get up to Python3 is as follows.

1. Check python version first just in case we can skip this

`which python3`


2. Install pyenv without using sudo

`curl https://pyenv.run | bash`


3. Add pyenv to shell init in ~/.bash.extensions

```
        export PYENV_ROOT="$HOME/.pyenv"
        export PATH="$PYENV_ROOT/bin:$PATH"
        eval "$(pyenv init -)"
        eval "$(pyenv virtualenv-init -)"

```
4. Source the ~/.bash.extensions

`source ~/.bash.extensions`


5. Install the new version of python and set it to the global version for your user:
```
      pyenv install 3.11.8
      pyenv global 3.11.8
```

6. Verify with `which` and `python --version`


## Install the scripts

There are three commands that I've found the most useful

`python3 -m pip install --user .`

`python3 -m pip install --user --upgrade .`

`python3 -m pip install --user --force-reinstall .`


Where the `.` must be the level of the script dir where the setup.py exists.
