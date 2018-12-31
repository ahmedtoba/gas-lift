from flask import Flask, render_template, request
import jinja2
from math import exp

app = Flask(__name__)
env = jinja2.Environment()
app.jinja_env.filters['zip'] = zip
env.globals.update(zip=zip)



#calculating Z_factor as a function of specific gravity, pressure and temperature
def Z_factor(g, p, t):
    T = t + 460
    Tpc = 168+325*g-12.5*g**2
    Ppc = 677+15*g-37.5*g**2

    Tr = T / Tpc
    Pr = p / Ppc
    t = Tpc / T


    a = -0.06125 * Pr * t * exp(-1.2*(1-t)**2)
    b = (14.76 * t - 9.76 * t**2 + 4.58 * t**3)
    c = 90.7 * t - 242.2 * t**2 + 42.4 * t**3
    d = 2.18 + 2.82 * t

    count = 1
    yi = 0.0125 * Pr * t * exp(-1.2*(1-t)**2)
    while True :
        y = yi
        f = a + (y + y**2 + y**3 + y**4)/(1-y)**3 - b * y**2 + c * y**d
        df = (1 + 4*y + 4*y**2 - 4*y**3 + y**4)/(1 - y)**4 - 2*b*y + c*d*y**(d-1)
        y = y - f/df
        check = abs(yi - y)
        if check < 10**-12:
            Y = abs(y)
            c = count
            break
        count += 1
        yi = y

    z = (0.06125 * Pr * t / Y )* exp(-1.2*(1-t)**2)
    return z


#Otis R table
portd = [0.125, 0.152, 0.173, 0.1875, 0.196, 0.2031, 0.213, 0.246, 0.25, 0.2812, 0.3125, 0.3437, 0.375, 0.4375, 0.5, 0.5625, 0.625, 0.75]
R_table = [[0.125, 0.1016, 0.0383, ''], [0.152, 0.1508, '', ''], [0.173, 0.1958, '', ''], [0.1875, '', 0.0863, 0.0359], [0.196, 0.2508, '', ''], [0.2031, '', 0.1013, ''], [0.213, 0.2966, '', ''], [0.246, 0.3958, '', ''], [0.25, '', 0.1534, 0.0638], [0.2812, '', 0.1942, ''], [0.3125, '', 0.2397, 0.0996], [0.3437, '', 0.29, ''], [0.375, '', 0.345, 0.1434], [0.4375, '', 0.4697, 0.1952], [0.5, '', '', 0.2562], [0.5625, '', '', 0.3227], [0.625, '', '', 0.3984], [0.75, '', '', 0.5738]]

def d_table(d):
    for x in portd:
        if x == d :
            return x
        elif round(x,4) == round(d,4):
            return x
        elif round(x,3) == round(d,3):
            return x
        elif round(x,2) == round(d,2):
            return x
        elif round(x,1) == round(d,1):
            return x

@app.route('/')

def index():
   result = False
   return render_template('index.html', result=result)

@app.route('/result',methods = ['POST', 'GET'])


def result():
   result = True
   if request.method == 'POST':
    #Extracting Values from html template forms
      Dp = float(request.form['Dp'])
      Di = float(request.form['Di'])
      Gs = float(request.form['gs'])
      gravity = float(request.form['y'])
      Tb = float(request.form['Tb'])
      Ts = float(request.form['Ts'])
      Pti = float(request.form['pti'])
      Pk = float(request.form['pk'])
      Pcs = float(request.form['ps'])
      Pwh = float(request.form['pwh'])
      ptm = float(request.form['ptm'])
      pcm = float(request.form['pcm'])
      qg = float(request.form['qg'])
      Pup = float(request.form['Pup'])
      Pdn = float(request.form['Pdn'])
      Tup = float(request.form['Tup'])
      k = float(request.form['k'])
      C = float(request.form['C'])
      S = float(request.form['S'])

    # calculating R or Ap/Ab of the valves:
      Ap = (qg)/(1248*C*Pup*(k/((k-1)*gravity*(Tup+460))*((Pdn/Pup)**(2/k)-(Pdn/Pup)**((k+1)/k)))**.5)
      dp = 1.1284*(Ap)**.5
      x = d_table(dp)
      row = portd.index(x)

      i = 1
      while i <= 3 :
          R = 0
          if R_table[row][i] != '' :
              R = R_table[row][i]
              break
          i += 1


      Phfd = Pwh + ptm
      Gfd = (Pti-Phfd)/Di
      Gt = (Tb-Ts)/Di

      valves = []
      d = (Pk-pcm-Pwh)/(Gs-(Pk-pcm)/40000)
      valves.append(d)
      while d < Di :
          d = (Pcs-pcm-Phfd+(Gs-Gfd)*d)/(Gs-(Pcs-pcm)/40000)
          valves.append(d)

      valves[-1] = Di
      Valves = [round(i) for i in valves]
      T = []
      DTP = []
      VOP = []
      DPD = []
      VCP = []
      DP60 = []
      TRO = []
      for l in Valves :
          T.append(round(Ts+Gt*l))
          DTP.append(round(Phfd+Gfd*l))
          VOP.append(round(Pcs*(1+l/40000)))

      for m, n  in zip(VOP,DTP):
          DPD.append(round((1-R)*m-S+R*n))

      for j, t in zip(DPD,T) :
          VCP.append(round(j+S*(1-R)))
          pd60 = round((520*j)/(t+460))
          Zs = Z_factor(gravity, pd60, 60)
          Zd = Z_factor(gravity, j, t)
          print(Zs,Zd)
          DP60.append(round((520*Zs*j)/((t+460)*Zd)))

      [TRO.append(round(d/(1-R)+S)) for d in DP60]


      return render_template("result.html", Valves=Valves, T=T, DTP=DTP, VOP=VOP, DPD=DPD, Pcs=round(Pcs), VCP=VCP, DP60=DP60, S=S, TRO=TRO, result=result, gs=Gs, y=gravity, Di=Di, Dp=Dp, pcm=pcm, ptm=ptm, pk=Pk, ps=Pcs, pwh=Pwh, Ts=Ts, Tb=Tb, pti=Pti, qg=qg, Pup=Pup, Pdn=Pdn, Tup=Tup, k=k, C=C )

if __name__ == '__main__':
   app.run(debug = True)
