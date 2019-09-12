% Main, ausschlie�lich Programmauswahl
%%
clear
clc
close all
%%
disp("Wirbelstroemungsprojekt, Numerik 1, Gruppe 118 \n")
pruef = true;

global Nx;
global sv;
global Re;
global obj;

while pruef == true
    Re = input("bei welcher Reynoldszahl sol die Str�mung simuliert werden? \n")
    Nx = input("wie gro� soll die Aufl�sung in X-Richtung sein?\n")
    sv = input("welches Seitenverh�ltnis soll vorliegen(x:y)?\n")
    a = input("Welche Str�mung soll dargestellt werden?\n 1) umstr�mtes Objekt\n 2) Lid driven cavitiy\n 3) Rohrstroemung\n 9) beenden\n")
    
    if a == 1
        umstroemung(sv,Nx,Re)
        pruef = false;
    elseif a == 2
        liddriven(sv,Nx,Re)
        pruef = false;
    elseif a == 9
        pruef = false;
    else
        disp("Antwort ung�ltig")
        pruef = true;
    end
end
