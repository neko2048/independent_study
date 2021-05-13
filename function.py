import sys
def progressbar(now, length, text):
        nbar = 10
        text = text + ' @ ' + str(now+1) + '/' + str(length)
        j = now / length # precent
        sys.stdout.write('\r')
        sys.stdout.write(f"[{'=' * int(nbar * j):{nbar}s}] {int(100 * j)}%  {text}")
        sys.stdout.flush()
        if now == length-1:
                sys.stdout.write('\r')
                sys.stdout.write(f"[{'=' * int(nbar * 1):{nbar}s}] {int(100 * 1)}%  {text}")
                print('\nend %s\n'%text)
