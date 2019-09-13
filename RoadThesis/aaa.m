clear all; clc; close all;

for (int i = 0; i < complexArray.Length / 2; i++) 
    outputArray[i] = 10.0 * Math.Log10((double)(Math.Sqrt((complexArray[i].Re * complexArray[i].Re) + (complexArray[i].Im * complexArray[i].Im))));