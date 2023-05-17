def boxymcboxface(message):
    """
    This function was obtained from `ARIBA GitHub`_

    It's just an aesthetics function which makes with a message inside. For example, with the word "message" it will write:
    
    .. code-block:: bash

       |==================================================|
       |                     example                      |
       |==================================================|


    :param message: Message to write inside the box of the function.
    :type messasge: string

    :return: an output text display box with a message inside.
    :rtype: output text display box.

    """
    str(message)
    #print('-' * 79)
    print ('\n')
    print('|', '=' * 50, '|', sep='')
    print('|', '{: ^48}'.format(message), '|')
    print('|', '=' * 50, '|', sep='')
    print ('\n')
    #print('-' * 79)