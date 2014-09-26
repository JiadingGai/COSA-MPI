close all;
clear all;

x = [32 16 8 4 3];
y = [4193 8389.44 17089 34721 51732];

for i=1:5
    %rs(i) = y(5) / y(i);
    rs(i) = y(6-i)
end

figure('Units','pixels', ...
        'Position', [100 100 500 375]);
h1 = plot(x,y);
%hold on;
%h2 = plot(x,rs);
set(h1,'Color',[0,0.5,0]);
set(h1                         , ...
  'Marker'          , 'o'         , ...
  'MarkerSize'      , 8           , ...
  'MarkerEdgeColor' , 'none'      , ...
  'MarkerFaceColor' , [.75 .75 1] );
set(h1,'LineWidth',2.0);

% set(h2,'Color',[0.5,0,0]);
% set(h2                         , ...
%   'Marker'          , 'o'         , ...
%   'MarkerSize'      , 8           , ...
%   'MarkerEdgeColor' , 'none'      , ...
%   'MarkerFaceColor' , [.75 .75 1] );
% set(h2,'LineWidth',2.0);


hTitle  = title ('Speed-up gain (parallel COSA)');
hXLabel = xlabel('Number of Processes (CPUs)'     );
hYLabel = ylabel('CPU Runtime (seconds)'          );
set([hTitle, hXLabel, hYLabel], ...
    'FontName'   , 'AvantGarde');
set([hXLabel, hYLabel]  , ...
    'FontSize'   , 10          );
set( gca                       , ...
    'FontName'   , 'Helvetica' );
set(gca             , ...
    'FontSize'   , 8           );
set( hTitle                    , ...
    'FontSize'   , 12          , ...
    'FontWeight' , 'bold'      );
set(gcf, 'PaperPositionMode', 'auto');
print -depsc2 finalPlot1.eps
print -dbmp finalPlot1.bmp
close;