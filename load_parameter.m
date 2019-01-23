function parameter = load_parameter(Settings)

% load 'Rasterraum informationen für Messung' 
start_position = regexp(Settings,'Rasterraum für Messung:');
end_position = regexp(Settings,'Schrittanzahl:');
parameter.rasterraum = regexp(Settings(start_position:end_position),'\d+(\.\d+)?','match');

% load 'Schrittanzahl' #Schritt = Schrittanzahl + 1
start_position = regexp(Settings,'Schrittanzahl:');
end_position = regexp(Settings,'Schrittweite');
parameter.schritt = regexp(Settings(start_position:end_position),'\d+(\.\d+)?','match');

% load 'Schrittweite'
start_position = regexp(Settings,'Schrittweite');
end_position = regexp(Settings,'Startposition Messung');
parameter.schrittweite = regexp(Settings(start_position:end_position),'\d+(\.\d+)?','match');

% load 'Koeffizienten'
start_position = regexp(Settings,'a5');
end_position = regexp(Settings,'a4');
koeffizienten= regexp(Settings(start_position:end_position),'[+-]?\d+(\.\d+)?','match');
parameter.koeffizienten(1,1) = koeffizienten(1,2);

start_position = regexp(Settings,'a4');
end_position = regexp(Settings,'a3');
koeffizienten= regexp(Settings(start_position:end_position),'[+-]?\d+(\.\d+)?','match');
parameter.koeffizienten(1,2) = koeffizienten(1,2);

start_position = regexp(Settings,'a3');
end_position = regexp(Settings,'a2');
koeffizienten= regexp(Settings(start_position:end_position),'[+-]?\d+(\.\d+)?','match');
parameter.koeffizienten(1,3) = koeffizienten(1,2);

start_position = regexp(Settings,'a2');
end_position = regexp(Settings,'a1');
koeffizienten= regexp(Settings(start_position:end_position),'[+-]?\d+(\.\d+)?','match');
parameter.koeffizienten(1,4) = koeffizienten(1,2);

start_position = regexp(Settings,'a1');
end_position = regexp(Settings,'a0');
koeffizienten= regexp(Settings(start_position:end_position),'[+-]?\d+(\.\d+)?','match');
parameter.koeffizienten(1,5) = koeffizienten(1,2);

start_position = regexp(Settings,'a0');
end_position = regexp(Settings,'Kalibrierte Temperatur');
koeffizienten= regexp(Settings(start_position:end_position),'[+-]?\d+(\.\d+)?','match');
parameter.koeffizienten(1,6) = koeffizienten(1,2);

% load 'Kalibrierte Temperatur'
start_position = regexp(Settings,'Kalibrierte Temperatur');
end_position = regexp(Settings,'Gaswerte');
parameter.temperatur = regexp(Settings(start_position:end_position),'\d+(\.\d+)?','match');

% load 'start_pos'
start_position = regexp(Settings,'Scanner-KOS');
end_position = regexp(Settings,'[Position');
parameter.start_pos = regexp(Settings(start_position:end_position),'[+-]?\d+(\.\d+)?','match');

% load 'sample_rate'
start_position = regexp(Settings,'Messeigenschaften');
end_position = length(Settings);
parameter.sample_rate = regexp(Settings(start_position:end_position),'\d+(\.\d+)?','match');

