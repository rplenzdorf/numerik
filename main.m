<<<<<<< HEAD
% Main, ausschlie�lich Programmauswahl

disp("Wirbelstroemungsprojekt, Numerik 1, Gruppe 118 \n")
pruef = true;

global Nx;
global sv;
global Re;
global obj;

while pruef == true
    Nx = input("wie gro� soll die Aufl�sung in X-Richtung sein?\n")
    sv = input("welches Seitenverh�ltnis soll vorliegen(x:y)?\n")
    a = input("Welche Str�mung soll dargestellt werden?\n 1) umstr�mtes Objekt\n 2) Lid driven cavitiy\n 9) beenden\n")
    
    if a==1
        b = input("welches Objekt soll umstr�mt werden?\n \n 1) Kreis\n 2) Rechteck\n 3) eigenes Objekt\n")
        if b == 1
            obj = (x-.5).^2+((y -.5)).^2<.2^2;
            wirbelstr(
        elseif b == 2
            obj = (x-.5).^10+((y -.5)).^10<.2^10;
        elseif b == 3
            obj = input("Objekt: \n")
        else
            disp("ung�ltige Eingabe")
        end
        pruef = false;
    elseif a == 2
        Re = input("bei welcher Reynoldszahl sol die Str�mung simuliert werden? \n")
        liddriven(sv,Nx,Re)
        pruef = false;
    elseif a == 9
        
    else
        disp("Antwort ung�ltig")
        pruef = true;
    end
end
=======
% Main, ausschlie�lich Programmauswahl
%%
clear
clc
close all
%%
fprintf("Wirbelstroemungsprojekt, Numerik 1, Gruppe 118\n\n")
pruef = true;

global Nx;
global sv;
global Re;
global obj;

while pruef == true
    Re = input("bei welcher Reynoldszahl sol die Str�mung simuliert werden? \n")
    Nx = input("wie gro� soll die Aufl�sung in X-Richtung sein?\n")
    sv = input("welches Seitenverh�ltnis soll vorliegen(x:y)?\n")
    a = input("Welche Str�mung soll dargestellt werden?\n 1) Umstr�mtes Objekt\n 2) Lid driven cavitiy\n 3) Rohrstroemung\n 4) Bewegtes Objekt\n 9) beenden\n")
    
    if a == 1
        umstroemung(sv,Nx,Re)
        pruef = false;
    elseif a == 2
        liddriven(sv,Nx,Re)
        pruef = false;
    elseif a == 3
        rohrstroemung(sv,Nx,Re)
        pruef = false;
    elseif a == 4
        bewegung(sv,Nx,Re)
        pruef = false;
    elseif a == 9
        pruef = false;
    else
        disp("Antwort ung�ltig")
        pruef = true;
    end
end
>>>>>>> a20f9a02db0b49fcf61f98de25989d73c04886c1
