@staticmethod
def QListWidget_to_list(qlist):
    """
    Method to convert QListWidget to python list

    Args:
        qlist (QListWidget)

    Returns:
        list: python list of strings
    """
    if qlist is not None:
        return [str(qlist.item(i).text()) for i in range(qlist.count())]
    else:
        return None
