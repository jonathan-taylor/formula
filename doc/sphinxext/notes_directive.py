import sys, os, cStringIO
try:
    from hashlib import md5
except ImportError:
    from md5 import md5

from docutils.parsers.rst import directives

import os, shutil, tempfile
import enthought.traits.api as traits

######################################################

NOTESPATH = "/home/jtaylo/bzr-repos/personal/courses/stats306B/notes"
WEBPATH = "http://stats306b.stanford.edu/restricted/"

class Notes(traits.HasTraits):

    pdf = traits.false
    name = traits.Str
    topic = traits.Str
    outfile = traits.Str

    def _name_changed(self):
        bdir = '%s/doc/www_new/restricted' % os.path.abspath(os.path.join(NOTESPATH, '..'))
        self.outfile = os.path.join("%s/%s.pdf" % (bdir, self.name))

    def __init__(self, name, topic, pdf=True):
        self.name = name
        self.topic = topic
        self.pdf = pdf

    def build_pdf(self):
        ndir = os.path.join(NOTESPATH, self.name)
        cmd = """
        cd %(ndir)s;
        pdflatex %(name)s_full.tex;
        rm %(name)s_full.aux;
        rm %(name)s_full.log;
        rm %(name)s_full.out;
        cp %(name)s_full.pdf %(outfile)s;
        """ % {'ndir':ndir, 'name':self.name,
               'outfile':self.outfile}
        print cmd
        os.system(cmd)

################################################

def notes_directive(name, arguments, options, content, lineno,
                    content_offset, block_text, state, state_machine):

    dirname = options['dirname']
    topic = options['topic']

    notes = []

    notes = Notes(dirname, topic)

    if not os.path.exists(notes.outfile) or 'force-latex' in options:
        notes.build_pdf()

    lines = ['Notes', '-----', 'Notes on %s available `here <restricted/%s.pdf>`_.' % (topic, dirname),'']
    
    state_machine.insert_input(
        lines, state_machine.input_lines.source(0))

    return []

options = {'dirname': directives.unchanged,
           'topic': directives.unchanged,
           'force-latex': directives.flag}
           

def setup(app):
    setup.app = app
    app.add_directive('notes', notes_directive, True, (0, 0, 0), **options)



