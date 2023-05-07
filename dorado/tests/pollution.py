def get_builtins():
    """Due to the way Python works, ``__builtins__`` can strangely be either a module or a dictionary,
    depending on whether the file is executed directly or as an import. I couldn’t care less about this
    detail, so here is a method that simply returns the namespace as a dictionary."""
    return getattr( __builtins__, '__dict__', __builtins__ )



G          = get_builtins()
G[ 'G' ]   = G
G[ 'axe' ] = axe

# maybe core instance should be elevated to builtin 