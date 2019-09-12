% Main, ausschließlich Programmauswahl
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
    Re = input("bei welcher Reynoldszahl sol die Strömung simuliert werden? \n")
    Nx = input("wie groß soll die Auflösung in X-Richtung sein?\n")
    sv = input("welches Seitenverhältnis soll vorliegen(x:y)?\n")
    a = input("Welche Strömung soll dargestellt werden?\n 1) umströmtes Objekt\n 2) Lid driven cavitiy\n 3) Rohrstroemung\n 9) beenden\n")
    
    if a == 1
        umstroemung(sv,Nx,Re)
        pruef = false;
    elseif a == 2
        liddriven(sv,Nx,Re)
        pruef = false;
    elseif a == 9
        pruef = false;
    else
        disp("Antwort ungültig")
        pruef = true;
    end
end
