function intersect(fig1, fig2)
  { 
  function step (a, b) {return a - a*Math.floor(a/b.length) + b.length*(Math.abs(a) - a)/2;};
  function rotate(A,B,C) {return (B.x-A.x)*(C.y-B.y)-(B.y-A.y)*(C.x-B.x)};
  function isCrossed(A,B,C,D){ return ((rotate(A,B,C)*rotate(A,B,D)) < 0 && (rotate(C,D,A)*rotate(C,D,B)) < 0);};
  function findCrP(A,B,C,D)
        {
          function findEq(E,F)
            {return {a: (F.y - E.y), b: (F.x - E.x), c: (F.x - E.x)*E.y - (F.y - E.y)*E.x};};
            diln = findEq(C,D).b*findEq(A,B).a - findEq(A,B).b*findEq(C,D).a;
            return {x: ((findEq(C,D).c*findEq(A,B).b - findEq(A,B).c*findEq(C,D).b)/diln), y: ((findEq(C,D).c*findEq(A,B).a - findEq(A,B).c*findEq(C,D).a)/diln)}
        };
  function isInPol(Pol, pt)
          {var c=0;
          for(var l=0; l<=Pol.length; l++)
            {
            if ((isCrossed(Pol[step(l, Pol)], Pol[step(l+1, Pol)], pt, {x: 500, y: pt.y + Math.pow(2,0.5)})))
                c++
            }
            return c%2;
          };
  function findS(mass)
    {
      for(var j = 0; j<mass.length-2; j++)
        {for(var n = j, Crp; n<mass.length-1; n++)
          {
            if (isCrossed(mass[j], mass[j+1], mass[n], mass[step(n+1, mass)]))
            {
              Crp = findCrP(mass[j], mass[j+1], mass[n], mass[step(n+1, mass)]);
              mass[j] = mass[n] = Crp;
            }
          }
        }
    for (var s=0, S=0; s<mass.length-1; s++)
    {
      S+=(mass[s].x*(mass[step(s+1, mass)].y-mass[step(s-1, mass)].y))/2 
    }
  return S}
  function AddCrossp(m1, m2)
    {
      function dist(A,B)
        {return Math.pow((Math.pow((A.x - B.x), 2) + Math.pow((A.y - B.y), 2)), 0.5)};
         var finFig = [];
        for (var i=0; i < m1.length; i++)
        { m1[i].visit = false;
          finFig.push(m1[i])
          var crossed = [];
          var k = 0; 
          for (var n=0; n < m2.length; n++)
           {m2[n].visit = false;
            if ((dist(m2[n], m1[i]) + dist(m2[n], m1[step(i+1, m1)])) == dist(m1[i], m1[step(i+1, m1)]))
              {if (isCrossed(m2[step(n-1, m2)], m2[step(n+1, m2)], m1[step(i, m1)], m1[step(i+1, m1)]))
                {crossed[k] = {x: m2[n].x, 
                      y: m2[n].y, 
                      alph: dist(m1[i], m2[n])/dist(m1[i],m1[step(i+1, m1)]), 
                      visit: false}
                k++;
                }
              else m2[n].EntOrEx = [0, 0];   
               }
            else 
              {
                if (isCrossed(m2[n], m2[step(n+1, m2)], m1[i], m1[step(i+1, m1)]))
                                    { var tx = findCrP(m2[n], m2[step(n+1, m2)], m1[i], m1[step(i+1, m1)]).x;
                                      var ty = findCrP(m2[step(n, m2)], m2[step(n+1, m2)], m1[i], m1[step(i+1, m1)]).y;
                                      crossed[k] = {x: tx, 
                                      y: ty, 
                                      alph:(dist({x: tx, y: ty}, m1[i])/dist(m1[i],m1[step(i+1, m1)])),
                                      visit: false
                                        };
                                      k++;
                                    }
              }
          }
          finFig = finFig.concat(crossed.sort(function(a,b) {if (a.alph > b.alph) return 1; if (a.alph <= b.alph) return -1}));
          }
      return finFig;
      };

      function marker(mass1, mass2)
      {
        function MarkTheWay(qpt, pt1, othm, qm, j)       
        {     
         if (("EntOrEx" in pt1) && (!((pt1.EntOrEx[0] == 0) && ((pt1.EntOrEx[1] == 0)))))
              {qpt.EntOrEx = [pt1.EntOrEx[1], pt1.EntOrEx[0]]; }
        else 
            {  var g = {x: (qpt.x + othm[step(j-1, othm)].x)/2, y: (qpt.y + othm[step(j-1, othm)].y)/2}
               qpt.EntOrEx = [isInPol(qm, g), (isInPol(qm, g) + 1)%2];
            }
         return qpt;
        };
        for (var i = 0; i < mass1.length; i++)
          {
            for (var n = 0; n < mass2.length; n++)
                  {if ((mass1[i].x == mass2[n].x) && (mass1[i].y == mass2[n].y))
                          { mass1[i].neig = n;
                            mass2[n].neig = i;
                            mass1[i] = MarkTheWay(mass1[i], mass1[step(i-1, mass1)], mass2, mass1, n);
                            mass2[n] = MarkTheWay(mass2[n], mass2[step(n-1, mass2)], mass1, mass2, i);
                          }
                  }
          }
        return [mass1, mass2];
        };
      function findRes(WalkMass, RefMass)
      {
        function walking(mass,d,t)
        {var Res = [];
          if (mass[t][d].EntOrEx[1] == 0)
            { d = d++;
              for(;(("EntOrEx" in mass[t][d])?(mass[t][d].EntOrEx[0] == 0):(isInPol(mass[(t+1)%2], mass[t][d]))); 
                d++) 
              { mass[t][d].visit = true;
                 Res.push(mass[t][d]);
                if ("neig" in mass[t][d]) {mass[(t+1)%2][mass[t][d].neig].visit = true;}
               }
              }
            else if (mass[t][d].EntOrEx[0] == 0)
            {d--
              for(;(("EntOrEx" in mass[t][d])?(mass[t][d].EntOrEx[1] == 0):(isInPol(mass[(t+1)%2], mass[t][d]))); 
                d--) 
              { mass[t][d].visit = true;
                Res.push(mass[t][d]);
                if ("neig" in mass[t][d]) {mass[(t+1)%2][mass[t][d].neig].visit = true;}
              }
            }
        return Res;
      };
    var Res = [];
    for (var q = 0, u, resulting = [], stPoint = {}; q < WalkMass.length; q++) 
      {
    if (WalkMass[q].visit == true)
      {continue;}
    else if (WalkMass[q].visit == false) 
      if ("EntOrEx" in WalkMass[q])
            { 
              stPoint = {x: WalkMass[q].x, y: WalkMass[q].y};
              WalkMass[q].visit = true;
              RefMass[WalkMass[q].neig].visit = true
              resulting = walking([WalkMass, RefMass], WalkMass[q].neig, 1);
              u = resulting[resulting.length - 1].neig;
              for (var i = 2, TempRes = []; 
                ((stPoint.x != resulting[resulting.length-1].x) ||
                 (stPoint.y != resulting[resulting.length-1].y)); 
                i++)
                { resulting = resulting.concat(walking([WalkMass, RefMass],u,(i%2)));
                  u = resulting[resulting.length - 1].neig;
                }
              if(findS(resulting)>0.001) {Res.push(resulting);}
          }
      else if (!("EntOrEx" in WalkMass[q]))
            {WalkMass[q].visit = true;}
    }  
    return [Res, RefMass, WalkMass];
    };
      if (rotate(fig1[0], fig1[1], fig1[2]) < 0)
        {fig1.reverse();}
      if (rotate(fig2[0], fig2[1], fig2[2]) < 0)
         {fig2.reverse();}
       var RESULT = [];
       var R = findRes(marker(AddCrossp(fig1, fig2), AddCrossp(fig2, fig1))[0],marker(AddCrossp(fig1, fig2), AddCrossp(fig2, fig1))[1]);
       var R1 = findRes(R[2], R[1]);
       RESULT = R[0].concat(R1[0]);
    return RESULT;
  };