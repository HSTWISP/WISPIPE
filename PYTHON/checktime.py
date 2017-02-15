import time


def checktime():
    now = time.strftime('%c')
    print 'Time check -  %s' % now



def main():
    checktime()


if __name__ == '__main__':
    main()
