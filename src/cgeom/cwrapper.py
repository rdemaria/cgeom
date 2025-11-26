import tinycwrap

clib = tinycwrap.CModule(__file__/"clib/base.c", __file__/"clib/path.c")