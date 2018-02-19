
C60 = ['D%s' % n for n in range(23, 27)]\
      + ['C%s' % n for n in range(22, 25)+range(27, 30)]\
      + ['C25T', 'C25B', 'C26T', 'C26B']\
      + ['B%s' % n for n in range(23, 28)]

def pack2sample(name):
    if name in C60: return 'C60'
    return 'Si'
