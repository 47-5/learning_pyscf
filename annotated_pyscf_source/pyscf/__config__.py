import os
import sys
import tempfile

#
# All parameters initialized before loading pyscf_conf.py will be overwritten
# by the dynamic importing procedure.
#

DEBUG = False

MAX_MEMORY = int(os.environ.get('PYSCF_MAX_MEMORY', 4000)) # MB
TMPDIR = os.environ.get('PYSCF_TMPDIR', tempfile.gettempdir())  # LRC NOTE: tempfile 是 Python 标准库里的一个模块，专门用来处理临时文件和临时目录。
                                                                # LRC NOTE: 这里的意思是：如果用户没有通过环境变量 PYSCF_TMPDIR 指定 PySCF 的临时目录，那就使用系统默认临时目录。
ARGPARSE = bool(os.getenv('PYSCF_ARGPARSE', False))

VERBOSE = 3  # default logger level (logger.NOTE)
UNIT = 'angstrom'

#
# Loading pyscf_conf.py and overwriting above parameters
#
for conf_file in (os.environ.get('PYSCF_CONFIG_FILE', None),  # LRC NOTE: os.environ是用于获取环境变量的函数
                  os.path.join(os.path.abspath('.'), '.pyscf_conf.py'),
                  os.path.join(os.environ.get('HOME', '.'), '.pyscf_conf.py')):
    if conf_file is not None and os.path.isfile(conf_file):
        break
else:
    conf_file = None

if conf_file is not None:
    with open(conf_file, 'r') as f:
        exec(f.read())
        """
        LRC Note: 
        f.read()返回类似 MAX_MEMORY = 8000\nVERBOSE = 4\nUNIT = 'angstrom'\nDUMPINPUT = False\ngto_mole_Mole_cart = True\n的东西，是一个长字符串
        exec()会把这个字符串当成python源码运行，这就相当于重新赋值。这就实现了覆盖默认值的行为
        """
    del f
del (os, sys, tempfile)

#
# All parameters initialized after loading pyscf_conf.py will be kept in the
# program.
#
