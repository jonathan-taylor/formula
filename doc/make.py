#!/usr/bin/env python
import fileinput
import glob
import os
import shutil
import sys

def check_build():
    build_dirs = ['build', 'build/doctrees', 'build/html', 'build/latex',
                  '_static', '_templates']
    for d in build_dirs:
        try:
            os.mkdir(d)
        except OSError:
            pass

def figs():
    os.system('cd pyplots/ && python make.py')

def html():
    check_build()
    if small_docs:
        options = "-D plot_formats=\"['png']\""
    else:
        options = ''
    if os.system('sphinx-build %s -P -b html -d build/doctrees . build/html' % options):
        raise SystemExit("Building HTML failed.")

    figures_dest_path = 'build/html/pyplots'
    if os.path.exists(figures_dest_path):
        shutil.rmtree(figures_dest_path)
#    shutil.copytree('pyplots', figures_dest_path)

def latex():
    check_build()
    #figs()
    if sys.platform != 'win32':
        # LaTeX format.
        if os.system('sphinx-build -b latex -d build/doctrees . build/latex'):
            raise SystemExit("Building LaTeX failed.")

        # Produce pdf.
        os.chdir('build/latex')

        # Copying the makefile produced by sphinx...
        if (os.system('pdflatex sampledoc.tex') or
            os.system('pdflatex sampledoc.tex') or
            os.system('makeindex -s python.ist Matplotlib.idx') or
            os.system('makeindex -s python.ist modMatplotlib.idx') or
            os.system('pdflatex sampledoc.tex')):
            raise SystemExit("Rendering LaTeX failed.")

        os.chdir('../..')
    else:
        print 'latex build has not been tested on windows'

def clean():
    for dirpath in ['build', 'examples']:
        if os.path.exists(dirpath):
            shutil.rmtree(dirpath)
    for pattern in ['mpl_examples/api/*.png',
                    'mpl_examples/pylab_examples/*.png',
                    'mpl_examples/pylab_examples/*.pdf',
                    'mpl_examples/units/*.png',
                    'pyplots/tex_demo.png',
                    '_static/matplotlibrc',
                    '_templates/gallery.html']:
        for filename in glob.glob(pattern):
            if os.path.exists(filename):
                os.remove(filename)

def all():
    #figs()
    html()
    latex()


funcd = {
    'figs'     : figs,
    'html'     : html,
    'latex'    : latex,
    'clean'    : clean,
    'all'      : all,
    }


small_docs = False

# Change directory to the one containing this file
current_dir = os.getcwd()
os.chdir(os.path.dirname(os.path.join(current_dir, __file__)))

if len(sys.argv)>1:
    if '--small' in sys.argv[1:]:
        small_docs = True
        sys.argv.remove('--small')
    for arg in sys.argv[1:]:
        func = funcd.get(arg)
        if func is None:
            raise SystemExit('Do not know how to handle %s; valid args are %s'%(
                    arg, funcd.keys()))
        func()
else:
    small_docs = False
    all()
os.chdir(current_dir)
os.system('killall -9 R')
