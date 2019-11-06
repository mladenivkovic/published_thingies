#!/usr/bin/env python3


#----------------------------------
# All stuff related to dealing
# with command line things. Be it
# reading arguments or interacting
#----------------------------------



#====================================
def yesno(*obj):
#====================================
    """
    Ask a yes/no question.
    Returns True if yes, False if no.

    You can use yesno() to print the question the same way
    you would use print() [without the formatting].
    """

    question = " ".join([str(o) for o in obj])

    while True:

        ans = input(question+' [y/n] ')
        try:
            form_ans = ans[0].lower()
        except IndexError:
            # happens for empty input
            form_ans = None

        if form_ans == 'y':
            return True
        elif form_ans == 'n':
            return False
        else:
            print("Didn't understand your input. Try again with 'y' or 'n'.")


    return







#===================================
def one_arg_present(default=None):
#===================================
    """
    Checks if there is exactly one cmd line arg present.
    If it is, it will be the return value.
    Otherwise, the 'default' parameter will be returned.
    """

    import sys

    if len(sys.argv) > 1:
        return sys.argv[1]
    else:
        return default
    
