# Import the .lvm Data : all files inside the select folder with .lvm will be loaded. 
(The lvm files inside the folder must arranged with measurement sequence)

# Load Parameters : instead input the parameters one by one, the parameter can also be loaded  with a .txt file (provided by Peter and Dennis)
(Important is the data sequence in the txt file:

---ALLGEMEIN---

G-Code Erstellungsdatum: 11-Jan-2019 ;
# Personen: Dennis Keller
# Zeitpunkt der Messung: 11:00
Dauer der Messung: 6 h ;

Bearbeitungskopf: Auﬂenmaﬂ X=290 ; Y=290 ; Innenmaﬂ X=260 ; Y=260 ;
# Spaltmaﬂ Bearbeitungskopf (Z): 0.6 mm ;


---KOORDINATEN---

Rasterraum f¸r Messung: 
X= ### mm ;
Y= ### mm ;
W= ### mm ;

Schrittanzahl:
X= ### ;
Y= ### ;
W= ### ;

Schrittweite (X,Y,W):
X= ### mm ;
Y= ### mm ;
W= ### mm ;

Startposition Messung (Maschinen-KOS):
X= ### mm ;
Y= ### mm ;
# Z= ### mm ; 
# W= ### mm ;

Startposition Messung (Scanner-KOS):
X= ### mm ; 
Y= ### mm ; 
# Z= ### mm ; [Position des Anemometer-Drahtes relativ zur Plattenoberseite]


---Anemometer---
Anemometerposition: 1 (X=299mm ; Y=299mm ; relativ zum Plattenursprung)
# Anemometer-Sonde: 15_15
Anemometer-Koeffizienten: a5=### ; a4=### ; a3=### ; a2=### ; a1=### ; a0=### ;

---Kalibrierte Temperatur---
Kalibrierte Temperatur: T=### ∞C;


---Gaswerte---
# Frischgas-Volumenstrom [Nl/min]: 25UK - 10BK - 5EEH ;
# Umw‰lz-Volumenstrom [m≥/h]: 36 ;
# Sauerstoff-Level [ppm]: 1000 ;

---Messeigenschaften---
# Messwerte pro Temperaturmesswert: ###;
)