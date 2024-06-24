#!/bin/bash

# Only shows the current directory on the prompt and not the default full path.
PS1='${debian_chroot:+($debian_chroot)}\[\033[01;32m\]\u@\h\[\033[00m\]:\[\033[01;34m\]\W\[\033[00m\]\$ '

# Shows directories in a different color than files.
alias ls="ls --color='auto'"
# Alias for a long and human readable form of ls.
alias ll='ls -lthG'

# Set the color scheme for ls.
LS_COLORS=$LS_COLORS:'di=0;35:' ; export LS_COLORS
