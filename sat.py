class Sat:
    t0 = 2451545.0

    def __init__(self, init=None, type='kep', t0=2451545.0, mu='geo'):
        self.t0 = t0
        if (mu == 'geo'):
            self.mu = 398600.4415  # км/с
        else:
            self.mu = 132712440018  # км/с

        if (type == 'kep'):
            if (init == None):
                init = [42164, 0, 0, 0, 0, 0]  # Нехай по умолчанию будет геостационарный спутник
            self.a, self.e, self.i, self.W, self.w, self.M0 = init
            self.cartesian()

        if (type == 'qv'):
            self.x, self.y, self.z, self.vx, self.vy, self.vz = init
            self.findOrbit()

    #  Процедуры вычисления:
    def findOrbit(self):
        q = [self.x, self.y, self.z]
        v = [self.vx, self.vy, self.vz]

        r = sqrt(dot(q, q))
        V = sqrt(dot(v, v))
        s = dot(q, v)
        h = cross(q, v)

        k = sqrt(self.mu)

        # Большая полуось
        a = 1 / abs(2. / r - V ** 2 / k ** 2)

        # Эксцентрисистет:
        e = sqrt( (s/k)**2 / a + (1 - r/a) ** 2)
        
        # Гиперболический случай:
        hyperbolic = True if e > 1 else False # Флаг гиперболического случая.
            

        if hyperbolic:
            dx = s / (k*sqrt(a))
            dy = (r/a + 1)
            H0 = arctanh(dx/dy)

        else:
            # Средняя аномалия:
            dy = s / (e * k * sqrt(a))
            dx = (a - r) / (a * e)
            E0 = arctan2(dy, dx)
            M0 = E0 - e * sin(E0)

        # Долгота восходящего узла:
        W = arctan2(h[0], -h[1])

        # Наклонение
        i = arctan2(sqrt(h[0] ** 2 + h[1] ** 2), h[2])

        # Аргумент перицентра:
        

        p = a*(e**2 - 1) if hyperbolic else a * (1 - e ** 2)

        dy = sqrt(p) * s
        dx = k * (p - r)
        vv = arctan2(dy, dx)

        if (sin(i) != 0):
            dy = self.z / sin(i)
            dx = self.x * cos(W) + self.y * sin(W)
            uu = arctan2(dy, dx)
        else:
            uu = 0

        w = uu - vv

        while (w < 0):
            w += 2 * pi
        if hyperbolic:
            self.a, self.e, self.i, self.W, self.w, self.H0 = a, e, i, W, w, H0
            return [a, e, i, W, w, H0]
        else:
            self.a, self.e, self.i, self.W, self.w, self.M0 = a, e, i, W, w, M0
            return [a, e, i, W, w, M0]

    def cartesian(self, t=t0, dt=None):
        a, e, i, W, w, M0 = self.get('kep')
        mu = self.mu
        t0 = self.t0
        # Поворотные матрицы:
        A = [[cos(W), -sin(W), 0],
             [sin(W), cos(W), 0],
             [0,      0,       1]]

        A = np.matrix(A)

        B = [[1, 0, 0], [0, cos(i), -sin(i)], [0, sin(i), cos(i)]]
        B = np.matrix(B)

        C = [[cos(w), -sin(w), 0], [sin(w), cos(w), 0], [0, 0, 1]]
        C = np.matrix(C)

        R = A * B * C  # Конечная поворотная матрица

        n = sqrt(mu * a ** (-3))

        if (dt):
            M = M0 + n*dt
        else:
            M = M0 + n * (t - t0)*86400

        E = M
        # Численное решение уравнения Кеплера
        while (abs(E - e * sin(E) - M) > 1e-12):
            E = E - (E - e * sin(E) - M) / (1 - e * cos(E))

        q = np.matrix([a * (cos(E) - e), a * sqrt(1 - e ** 2) * sin(E), 0])
        dq = np.matrix([-a * n * sin(E) / (1 - e * cos(E)), a * sqrt(1 - e ** 2) * cos(E) * n / (1 - e * cos(E)), 0])

        q = R * q.T  # T - преобразование строки к столбцу, суть транспозиция
        dq = R * dq.T

        self.x, self.y, self.z = q.A1  # A1 - преобразование к человеческому типу
        self.vx, self.vy, self.vz = dq.A1
        return [self.x, self.y, self.z, self.vx, self.vy, self.vz]

    # Процедуры установки новых параметров
    def set(self, settings, type='qv'):
        if (type == 'qv'):
            self.x, self.y, self.z, self.vx, self.vy, self.vz = settings
            self.findOrbit()
        if (type == 'kep'):
            self.a, self.e, self.i, self.W, self.w, self.M0 = settings
            self.cartesian()

    # Процедуры возвращающие параметры:
    def get(self, type='qv'):
        if (type == 'qv'):
            return [self.x, self.y, self.z, self.vx, self.vy, self.vz]

        if (type == 'kep'):
            return [self.a, self.e, self.i, self.W, self.w, self.M0]

    def subsat(self, JD=None, dt=None):
        if (JD == None):
            JD = self.t0 + dt/86400

        H = siderial(JD)  # Вычисляем звёздное время и матрицу поворота A

        A = np.matrix([[cos(H),     sin(H),     0],
                       [-sin(H),    cos(H),     0],
                       [0,          0,          1]])

        q = self.cartesian(JD)
        x = np.matrix(q[0:3])

        y = np.array(A*x.T)
        y_norm = sqrt(y[0]**2 + y[1]**2 + y[2]**2)

        L = float(arctan2(y[1], y[0])) # float для того, чтобы не возвращался объект типа matrix
        phi = float(arcsin(y[2]/y_norm))

        return [L, phi]

