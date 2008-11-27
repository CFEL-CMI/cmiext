#!/bin/tcsh

python setup.py sdist --dist-dir=new_distribution --formats=gztar \
    && scp new_distribution/jkext*.tar.gz ssh.rz-berlin.mpg.de:public_html/uploads/Computer/jkext-latest.tar.gz \
    && ssh ssh.rz-berlin.mpg.de chmod 644 /home/jochen/public_html/uploads/Computer/jkext-latest.tar.gz \
    && rm new_distribution/jkext*.tar.gz \
    && rmdir new_distribution
