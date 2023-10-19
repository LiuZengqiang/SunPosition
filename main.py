#!/usr/bin/env python3
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys

if __name__ == '__main__':
    filename = "./input.csv"
    if (len(sys.argv) > 1):
        filename = str(sys.argv[1])

    # 高度角为与地面平面的夹角
    # 方位角为以正北为0度角
    fig = plt.figure(1)
    ax = fig.add_axes(Axes3D(fig))

    input_file = open(filename)

    nrows = 0
    station = []
    year = []
    month = []
    day = []
    hour = []
    min = []
    sec = []
    lon = []
    lat = []
    TimeZone = []

    for idx, _ in enumerate(input_file.readlines()):
        if (idx == 0):
            pass
        else:
            temp = _.rstrip().split(',')
            nrows += 1
            station.append(temp[0])
            year.append(int(temp[1]))
            month.append(int(temp[2]))
            day.append(int(temp[3]))
            hour.append(int(temp[4]))
            min.append(int(temp[5]))
            sec.append(int(temp[6]))
            lon.append(float(temp[7]))
            lat.append(float(temp[8]))
            TimeZone.append(float(temp[9]))

    altitude = np.array([])
    azimuth = np.array([])
    x = np.array([])
    y = np.array([])
    z = np.array([])

    for n in range(1, nrows):
        m = n - 1
        # 年积日的计算
        # 儒略日 Julian day(由通用时转换到儒略日)
        JD0 = int(365.25 * (year[m] - 1)) + int(30.6001 *
                                                (1 + 13)) + 1 + hour[m] / 24 + 1720981.5

        if month[m] <= 2:
            JD2 = int(365.25 * (year[m] - 1)) + int(30.6001 *
                                                    (month[m] + 13)) + day[m] + hour[m] / 24 + 1720981.5
        else:
            JD2 = int(365.25 * year[m]) + int(30.6001 *
                                              (month[m] + 1)) + day[m] + hour[m] / 24 + 1720981.5

        # 年积日 Day of year
        DOY = JD2 - JD0 + 1

        # N0   sitar=θ
        N0 = 79.6764 + 0.2422 * (year[m] - 1985) - int((year[m] - 1985) / 4.0)
        sitar = 2 * math.pi * (DOY - N0) / 365.2422
        ED1 = 0.3723 + 23.2567 * math.sin(sitar) + 0.1149 * math.sin(2 * sitar) - 0.1712 * math.sin(
            3 * sitar) - 0.758 * math.cos(sitar) + 0.3656 * math.cos(2 * sitar) + 0.0201 * math.cos(3 * sitar)
        ED = ED1 * math.pi / 180  # ED本身有符号

        if lon[m] >= 0:
            if TimeZone == -13:
                dLon = lon[m] - (math.floor((lon[m] * 10 - 75) / 150) + 1) * 15.0
            else:
                dLon = lon[m] - TimeZone[m] * 15.0  # 地球上某一点与其所在时区中心的经度差
        else:
            if TimeZone[m] == -13:
                dLon = (math.floor((lon[m] * 10 - 75) / 150) + 1) * 15.0 - lon[m]
            else:
                dLon = TimeZone[m] * 15.0 - lon[m]
        # 时差
        Et = 0.0028 - 1.9857 * math.sin(sitar) + 9.9059 * math.sin(2 * sitar) - 7.0924 * math.cos(
            sitar) - 0.6882 * math.cos(2 * sitar)
        gtdt1 = hour[m] + min[m] / 60.0 + sec[m] / 3600.0 + dLon / 15  # 地方时
        gtdt = gtdt1 + Et / 60.0
        dTimeAngle1 = 15.0 * (gtdt - 12)
        dTimeAngle = dTimeAngle1 * math.pi / 180
        latitudeArc = lat[m] * math.pi / 180

        # 高度角计算公式
        HeightAngleArc = math.asin(
            math.sin(latitudeArc) * math.sin(ED) + math.cos(latitudeArc) * math.cos(ED) * math.cos(dTimeAngle))
        # 方位角计算公式
        CosAzimuthAngle = (math.sin(HeightAngleArc) * math.sin(latitudeArc) - math.sin(ED)) / math.cos(
            HeightAngleArc) / math.cos(latitudeArc)
        AzimuthAngleArc = math.acos(CosAzimuthAngle)
        HeightAngle = HeightAngleArc * 180 / math.pi
        ZenithAngle = 90 - HeightAngle
        AzimuthAngle1 = AzimuthAngleArc * 180 / math.pi

        if dTimeAngle < 0:
            AzimuthAngle = 180 - AzimuthAngle1
        else:
            AzimuthAngle = 180 + AzimuthAngle1
        altitude = np.append(altitude, HeightAngle)
        azimuth = np.append(azimuth, AzimuthAngle)
        aa = np.cos(HeightAngle / 180 * np.pi) * np.cos(AzimuthAngle / 180 * np.pi)
        cc = np.sin(HeightAngle / 180 * np.pi)
        bb = np.cos(HeightAngle / 180 * np.pi) * np.sin(AzimuthAngle / 180 * np.pi)
        ll = aa * aa + bb * bb + cc * cc
        ll = np.sqrt(ll)
        aa = aa / ll
        bb = bb / ll
        cc = cc / ll
        x = np.append(x, aa)
        y = np.append(y, bb)
        z = np.append(z, cc)
        if (min[m] % 5 == 0):
            print('站位：' + station[m] + ' 时间：%d:%d 高度角(deg): %f 方位角(deg): %f ' %
                  (hour[m], min[m], HeightAngle, AzimuthAngle))

    ax.plot(x, y, z)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.show()
