# -*- coding: utf-8 -*-
"""
Created on Fri Jun 17 13:19:16 2022

@author: DELL
"""

# =============================================================================
# #format of class method
#
# class Foo:
# 
#     def __init__(self, name):
#         self.name = name
# 
#     def ord_func(self):
#         """ 定义普通方法，至少有一个self参数 """
# 
#         # print self.name
#         print('普通方法' )
# 
#     @classmethod
#     def class_func(cls):
#         """ 定义类方法，至少有一个cls参数 """
# 
#         print('类方法')
#         
#     @staticmethod
#     def static_func():
#         """ 定义静态方法 ，无默认参数"""
# 
#         print('静态方法') 
# 
# 
# #方法的定义和使用
# # 调用普通方法
# Foo.ord_func()
# 
# # 调用类方法
# Foo.class_func()
# 
# # 调用静态方法
# Foo.static_func()
# =============================================================================
print('/////////////example 1 ://///////////')
# example 1
# 普通方法
class Student:
    def __init__(self, name, score1, score2):
        self.name = name
        self.score1 = score1
        self.score2 = score2
    def add_score(self):
        return self.score1 + self.score2#注
student = Student("lihua", 99, 95)#实例化
sname = student.name #>>"lihua"
student.score1
student.score2
sums = student.add_score()

print(sname)
print(sums)#>>194

print('/////////////example 2 ://///////////')
# example 2
# class方法
class Student:
    school = '田村小学'
    @classmethod
    def add_score(cls,name,score1,score2):
        print('%s的总成绩是 %s ' %(name, score1 + score2))
        print('%s的学校是 %s' %(name, cls.school))
Student.add_score('李华',91,80)
#>>李华的总成绩是 170 
#>>李华的学校是 田村小学
print('/////////////example 3 ://///////////')

# example 3
# 静态方法
# 在def前面加上**@staticmethod**，这种类方法是静态的类方法，他的一个特点是参数可以为空，
# 通过类名去调用，不需要self，就将这个方法当成一个普通的函数使用。
class Student:
    @staticmethod
    def add_score(school,name,score1,score2):
        print('%s的总成绩是 %s ' %(name, score1 + score2))
        print('%s的学校是 %s' %(name, school))
Student.add_score('汇师小学','Li',81,80)
#>>李华的总成绩是 170 
#>>李华的学校是 田村小学
