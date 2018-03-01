import numpy as np

from numpy import sin, cos, sqrt, arctan2, arcsin
from .consts import *
import re

def reduce(angle):
	"""
	Приводит угол к интервалу [0, 2pi)
	:param angle:
	:return:
	"""
	while (angle < 0):
		angle += 2*pi

	while (angle >= 2*pi):
		angle -= 2*pi

	return  angle


def siderial(JD):
	"""
	Вычисление звёздного времени на некоторый момент.
	:param JD: Юлианская дата момента
	:return:
	"""
	Tu = (JD - JD2000)/36525
	t = (JD % 1) - 0.5

	H0 = 24110.54841 + 8640184.812866 * Tu + 0.093104 * Tu**2 - 6.21e-6 * Tu**3
	s = (H0 / 86400 + t) * 2*pi
	# Приведение к первому периоду:


	return s


def readeph(str, strtype, fileformat='MEGNO'):
	"""
	Парсинг выходных файлов численной модели движения ИСЗ
	:param str: Строка из файла
	:param strtype: Тип строки (q - коорд, v - скорости, date - строка с датой)
	:param fileformat: Формат .EPH файла - с мегно или без.
	:return: Массив значений или None в случае неудачи
	"""

	if fileformat == 'MEGNO':
		# ------------#-----x--------------------y--------------------z--------------------MEGNO-----------------
		qrx = r'\s+(\d+)\s+([+-]?\d+\.\d+E?[+-]?\d*)\s+([+-]?\d+\.\d+E?[+-]?\d*)\s+([+-]?\d+\.\d+E?[+-]?\d*)\s+([+-]?\d+\.\d+E?[+-]?\d*)\s*'

		# ----------vx-------------------vy-------------------vz------------------ mMEGNO-------------
		vrx = r'\s+([+-]?\d+\.\d+E?[+-]?\d*)\s+([+-]?\d+\.\d+E?[+-]?\d*)\s+([+-]?\d+\.\d+E?[+-]?\d*)\s+([+-]?\d+\.\d+E?[+-]?\d*)\s*'

		# --------------JD0---------dt---------------yr------mnth-----dy------HR------MIN------SEC--------
		daterx = r'\s?(\d+\.\d?)\s+(\d+\.\d+)\s+\(\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)H\s+(\d+)M\s+(\d+\.\d+)S\)\s?'

	else:
		# ------------#-----x--------------------y--------------------z--------------------MEGNO-----------------
		qrx = r'\s*([+-]?\d+\.\d+E?[+-]?\d*)\s+([+-]?\d+\.\d+E?[+-]?\d*)\s+([+-]?\d+\.\d+E?[+-]?\d*)\s*'
		# ----------vx-------------------vy-------------------vz------------------ mMEGNO-------------
		vrx = r'\s*([+-]?\d+\.\d+E?[+-]?\d*)\s+([+-]?\d+\.\d+E?[+-]?\d*)\s+([+-]?\d+\.\d+E?[+-]?\d*)\s*'
		# --------------JD0---------dt---------------yr------mnth-----dy------HR------MIN------SEC--------
		daterx = r'\s?(\d+\.\d?)\s+(\d+.\d+)\s+TT \(\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)H\s+(\d+)M\s+(\d+\.\d+)S UTC\)\s?'

	if (strtype == 'date'):
		match = re.match(daterx, str)
	elif (strtype == 'q'):
		match = re.match(qrx, str)
	elif (strtype == 'v'):
		match = re.match(vrx, str)
	else:
		return None

	if (match):  # Проверяем, есть ли совпадение:
		out = []  # Создаём пустой массив чтобы записать в него float и вернуть
		for s in match.groups():
			out.append(float(s))  # Преобразуем из str в float
		return out
	else:
		return None


def subsat(sat, JD=None, dt=0):
	"""
	Поиск подспутниковой точки (только для спутников Земли)
	:param sat: Объект типа Sat
	:param JD: Юлианская дата на момент прогноза
	:param dt: ИЛИ смещение относительно sat.t0
	:return: list(широта, долгота) в радианах
	"""
	if (JD == None):
		JD = sat.t0 + dt / 86400

	H = siderial(JD)  # Вычисляем звёздное время и матрицу поворота A

	A = np.matrix([[cos(H),  sin(H),	0],
				   [-sin(H), cos(H),	0],
				   [	  0,	  0,	1]])

	q = sat.cartesian(JD)
	x = np.matrix(q[0:3])

	y = np.array(A * x.T)
	y_norm = sqrt(y[0] ** 2 + y[1] ** 2 + y[2] ** 2)

	L = float(arctan2(y[1], y[0]))  # float для того, чтобы не возвращался объект типа matrix
	phi = float(arcsin(y[2] / y_norm))

	return [L, phi]


def ephem(JD, target):
	# * Получение координат и скорстей в удобоваримом виде и единицах
	# * Пересчёт координат относительно Солнца, а не барицентра СС
	# * Вычисление координат Земли, а не барицентра ЗЛ

	qv = eph.position_and_velocity('sun', JD)
	qsun = qv[0]
	vsun = qv[1]

	target = target.lower() # Чтобы не промахнуться с регистром

	if (target == 'earth'): # Достаём координаты Земли из барицентра
		qv = eph.position_and_velocity('moon', JD)
		qmoon = qv[0]
		vmoon = qv[1]

		qv = eph.position_and_velocity('earthmoon', JD)
		q = qv[0] - (1./(1 + EMRAT))*qmoon
		v = qv[1] - (1./(1 + EMRAT))*vmoon

	else:
		qv = eph.position_and_velocity(target, JD)
		q = qv[0]
		v = qv[1]

	if (target != 'earth'):  # Луна в геоцентрической системе
		q -= qsun
		v -= vsun

	v /= 86400 # Приводим к км/с

	# Магия работы с эфемеридами:
	q = q.transpose()[0]
	v = v.transpose()[0]
	return np.array([q[0], q[1], q[2], v[0], v[1], v[2]])