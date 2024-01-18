nly shows the current directory on the prompt and not the default full path.
PS1='${debian_chroot:+($debian_chroot)}\[\033[01;32m\]\u@\h\[\033[00m\]:\[\033[01;34m\]\W\[\033[00m\]\$ '

# Alias for a long and human readable form of ls.
alias ll='ls -lthG'

# Shows directories in a different color than files.
LS_COLORS=$LS_COLORS:'di=0;35:' ; export LS_COLORS


