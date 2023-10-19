## 用于计算太阳位置的脚本  

### 简介  
一个用于根据`地点`、`时间`计算太阳位置的脚本
>> 计算结果可以跟[SunCalc sun position](https://www.suncalc.org)进行对比验证

### 使用  
#### 输入  
``*.csv``文件，里面包含 ``station, year, month, day, hour, min, sec, lon, lat, TimeZone``共10列数据，表示待求的``地点, 年, 月, 日, 小时, 分钟, 秒, 经度, 纬度, 时区``。   

#### 输出  
打印出待求地点、时间的太阳高度角、方位角（角度制）.  
高度角：``地点-太阳``向量与地面的夹角大小，[0,90]，平行地面为0，垂直地面为90.
方位角：``地点-太阳``向量与正北方的夹角大小，[0,360]，指向正北为0(360)，正东为90，正南为180，正西为270.

#### 示例  
默认以``input.csv``作为输入
```shell
main.py
```
运算结果(部分)：
```
站位：北京 时间：6:0 高度角(deg): -4.521586 方位角(deg): 86.663832 
站位：北京 时间：6:5 高度角(deg): -3.563963 方位角(deg): 87.469214 
站位：北京 时间：6:10 高度角(deg): -2.605746 方位角(deg): 88.273048 
站位：北京 时间：6:15 高度角(deg): -1.647124 方位角(deg): 89.075786 
站位：北京 时间：6:20 高度角(deg): -0.688285 方位角(deg): 89.877879 
...
```