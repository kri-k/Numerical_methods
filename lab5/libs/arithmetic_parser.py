import re


class FunctionFromStr:
    # use exec and regexp coz it's easier =)
    # may be should convert to postfix notation
    def __init__(self, s):
        self.valid = FunctionFromStr.validate(s)
        self.s_func = s if self.valid else ''
        self.args = []

    @staticmethod
    def validate(s):
        # available_funcs = '|'.join(['sin', 'cos', 'exp'])
        # available_ops = '|'.join(['\+', '-', '\*', '\/', '\^'])
        # func_pat = '(({0})\(\w+\))'.format(available_funcs)
        # func_var_const_pat = '({0}|\w+)'.format(func_pat)
        # pat = '^({0}\s*{1}\s*)*{0}$'.format(func_var_const_pat, available_ops)
        # pat = re.compile(pat)
        # return pat.match(s) is not None
        return True

    def set_args_name(self, *args):
        self.args = args

    def __call__(self, *args, **kwargs):
        s = '''
from math import *
answer = ({})'''.format(self.s_func)
        for i in range(min(len(self.args), len(args))):
            kwargs[self.args[i]] = args[i]
        try:
            exec(s, None, kwargs)
        except (SyntaxError, NameError) as e:
            self.valid = False
            return None
        return kwargs['answer']
