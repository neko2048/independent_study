import sys
import numpy as np
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

def buoyancy(th, qv):
        """
        input: th(z, y, x) [K], qv(z, y, x) [kg/kg]
        retrun buoyancy(z, y, x) [m/s^2]
        """
        thv = th * (1 + 0.622 * qv) # theta_v
        thvbar = np.mean(thv, (1, 2))
        thvbar = np.tile(thvbar[:, np.newaxis, np.newaxis], (1, 128, 128)) # horizontal averaged theta_v
        B = (thv - thvbar) / thvbar
        return B
