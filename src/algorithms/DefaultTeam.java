package algorithms;

import java.awt.Point;
import java.util.ArrayList;

public class DefaultTeam {

	
	public ArrayList<Point> calculAngularTSP(ArrayList<Point> points, int edgeThreshold, ArrayList<Point> hitPoints) {
	    int[][] paths=new int[points.size()][points.size()];
	    double[][] dist=new double[points.size()][points.size()];
	    
	    floydWarshall(points, edgeThreshold, paths, dist);
	    ArrayList<Point> result = calculATSP(points, edgeThreshold, hitPoints, paths, dist);
	    
	    
	    
	    ArrayList<Point> tmp;
	    for(int i = 0; i < 50; i++) {
	    	tmp = calculATSP(points, edgeThreshold, hitPoints, paths, dist);
	    	if(Evaluator.score(tmp) < Evaluator.score(result)) {
	    		result = tmp;
	    	}
	    }
	    
	    //return result;
	    
	    System.out.println("SCORE NORMAL: "+Evaluator.score(result));
	    while (Evaluator.score(result)>Evaluator.score(bruteforceWindow(result, edgeThreshold, dist))) result=bruteforceWindow(result, edgeThreshold, dist);
	    System.out.println("SCORE BRUTEFORCE: "+Evaluator.score(result));
	    while (Evaluator.score(result)>Evaluator.score(localSearchCross(result, edgeThreshold, paths, dist, points))) result=localSearchCross(result, edgeThreshold, paths, dist, points);
	    System.out.println("SCORE LOCAL CROSS: "+Evaluator.score(result));
	    return result;
  }
  
  public static ArrayList<ArrayList<Point>> permutation(ArrayList<Point> topermut){
	  return tmppermutation(topermut, new ArrayList<>(), new ArrayList<>());
  }
  
  
  
  public static ArrayList<ArrayList<Point>> tmppermutation(ArrayList<Point> topermut, ArrayList<Point> permuting, ArrayList<ArrayList<Point>> list){
	  if(topermut.isEmpty()) {
		  list.add(permuting);
	  }
	  ArrayList<Point> tmptopermut, tmppermuting;
	  
	  for(int i = 0; i < topermut.size(); i++) {
		  tmptopermut = new ArrayList<>(topermut);
		  tmppermuting = new ArrayList<>(permuting);
		  tmppermuting.add(tmptopermut.remove(i));
		  tmppermutation(tmptopermut, tmppermuting, list);
	  }
	  
	  return list;
  }
  
  
  public void validpermutation(ArrayList<ArrayList<Point>> list, Point p, Point q, int edgeThreshold, double[][] dist) {
	  for(ArrayList<Point> tmp: list) {
		  tmp.add(0, p);
		  tmp.add(q);
	  }
	  
	  ArrayList<ArrayList<Point>> notvalid = new ArrayList<>();
	  ArrayList<Point> tmp;
	  for(int i = 0; i < list.size(); i++) {
		  tmp = list.get(i);
		  for(int j = 0; j < tmp.size()-1; j++) {
			  if(dist[j][j+1]> edgeThreshold) {
				  notvalid.add(tmp);
				  break;
			  }
		  }
	  }
	  
	  list.removeAll(notvalid);
	  
  }
  
  public ArrayList<Point> minimumcost(ArrayList<ArrayList<Point>> list){
	  ArrayList<Point> min = null;
	  double minscore = Double.MAX_VALUE;
	  double tmpscore;
	  for(ArrayList<Point> tmp: list) {
		  if(min == null) {
			  min = list.get(0);
			  minscore = minscore(min);
		  }
		  tmpscore = minscore(tmp);
		  if(tmpscore < minscore) {
			  min = tmp;
		  }
	  }
	  return min;
	  
  }

  public double minscore(ArrayList<Point> points){
	  double s=0;
      for (int i=0;i<points.size()-1;i++) s+=points.get(i).distance(points.get((i+1)%points.size()));
      
	  double scal;
	  Point p,q,r;
	  FloatPoint v1, v2;
	  double value;
	  for(int i = 1; i <= points.size()-1; i++) {
		  p = points.get((i-1)%points.size());
		  q = points.get((i)%points.size());
		  r = points.get((i+1)%points.size());
		  v1 = new FloatPoint(q.getX()-p.getX(),q.getY()-p.getY());
		  v2 = new FloatPoint(r.getX()-q.getX(),r.getY()-q.getY());
		  scal = v1.x*v2.x + v1.y*v2.y;
		  value = scal/(v1.dist()*v2.dist());
		  if(value < -1) {
			  value = -1;
		  }else if(value > 1) {
			  value = 1;
		  }
		  s += Evaluator.angle(p, q, r);
	  }
	  return s;
  }
  
  public ArrayList<Point> bruteforceWindow(ArrayList<Point> points, int edgeThreshold, double[][] dist){
	  Point p1, p2;
	  Point i1, i2, i3, i4;
	  ArrayList<ArrayList<Point>> permut;
	  ArrayList<Point> topermut;
	  ArrayList<Point> sol;
	  ArrayList<Point> res;
      for (int i=0;i<points.size();i++){
    	  sol = null;
    	  p1 = points.get(Math.abs((i-1)%points.size()));
    	  i1 = points.get(i);
    	  i2 = points.get((i+1)%points.size());
    	  i3 = points.get((i+2)%points.size());
    	  i4 = points.get((i+3)%points.size());
    	  p2 = points.get((i+4)%points.size());
    	  topermut = new ArrayList<>();
    	  topermut.add(i1);
    	  topermut.add(i2);
    	  topermut.add(i3);
    	  topermut.add(i4);
    	  
    	  permut = permutation(topermut);
    	  validpermutation(permut, p1, p2, edgeThreshold, dist);
    	  sol = minimumcost(permut);

    	  if(sol != null) {
    		  res = new ArrayList<>();
    		  res.addAll(sol);
    		  for(int k = 0; k < points.size()-6; k++) {
    			  res.add(points.get((i+5+k)%points.size()));
    		  }
    		  return res;
    	  }
      }
      return points;
  }
  
  public ArrayList<Point> calculATSP(ArrayList<Point> points, int edgeThreshold, ArrayList<Point> hitPoints, int[][] paths, double[][] dist){
	   ArrayList<Point> result = new ArrayList<Point>();
	   ArrayList<Point> dup = new ArrayList<>(hitPoints);
	    
	    
	    
	    Point m, n = null;
	    int select = (int)(Math.random()*dup.size());
	    select = 0;
	    m = dup.remove(select);
	    while(dup.size() > 0) {
	    	n = dup.get(0);
	    	for(int j = 1; j < dup.size(); j++) {
	    		if(dist[points.indexOf(m)][points.indexOf(n)] > dist[points.indexOf(m)][points.indexOf(dup.get(j))]) {
	    			n = dup.get(j);
	    		}
	    	}
	    	result.addAll(expandPaths(points, m, n, paths, dist));
	    	m=dup.remove(dup.indexOf(n));
	    }
	    result.addAll(expandPaths(points, n, result.get(0), paths, dist));
	    
	    return result;
  }
  
  
  
  
  public void floydWarshall(ArrayList<Point> points, int edgeThreshold, int[][] paths, double[][] dist) {
	    for (int i=0;i<paths.length;i++) for (int j=0;j<paths.length;j++) paths[i][j]=i;


	    for (int i=0;i<paths.length;i++) {
	      for (int j=0;j<paths.length;j++) {
	        if (i==j) {dist[i][i]=0; continue;}
	        if (points.get(i).distance(points.get(j))<=edgeThreshold) dist[i][j]=points.get(i).distance(points.get(j)); else dist[i][j]=Double.POSITIVE_INFINITY;
	        paths[i][j]=j;
	      }
	    }

	    for (int k=0;k<paths.length;k++) {
	      for (int i=0;i<paths.length;i++) {
	        for (int j=0;j<paths.length;j++) {
	          if (dist[i][j]>dist[i][k] + dist[k][j]){
	            dist[i][j]=dist[i][k] + dist[k][j];
	            paths[i][j]=paths[i][k];

	          }
	        }
	      }
	    }
  }
  
  //ajoute pas le dernier
  public ArrayList<Point> expandPaths(ArrayList<Point> points, Point p, Point q, int[][] paths, double[][] dist){
	  ArrayList<Point> res = new ArrayList<>();
	  int i = points.indexOf(p);
	  int j = points.indexOf(q);
	  
	  while(i!=j) {
	    	p=points.get(i);       

	        res.add(p);

	        i=paths[i][j];
	        
	 }
	 
	  
	  return res;
  }
  

  

  
  private ArrayList<Point> localSearchCross(ArrayList<Point> points, int edgeThreshold, int[][] paths, double[][] dist, ArrayList<Point> original){
      for (int i=0;i<points.size();i++){
          for (int j=i+2;j<points.size() ;j++){
        	  double a=dist[i][(i+1)%points.size()];
              double b=dist[j%points.size()][(j+1)%points.size()];
              double c=dist[i][j%points.size()];
              double d=dist[(i+1)%points.size()][(j+1)%points.size()];
              
              double a1;
              double a2;
              
              a1 = Evaluator.angle(points.get((i-1+points.size())%points.size()), points.get((i)%points.size()), points.get((i+1)%points.size())) +
            		  Evaluator.angle(points.get((i)%points.size()), points.get((i+1)%points.size()), points.get((i+2)%points.size())) +
            		  Evaluator.angle(points.get((j-1+points.size())%points.size()), points.get((j)%points.size()), points.get((j+1)%points.size())) +
            		  Evaluator.angle(points.get((j)%points.size()), points.get((j+1)%points.size()), points.get((j+2)%points.size()));
              
              a2 = Evaluator.angle(points.get((i-1+points.size())%points.size()), points.get((i)%points.size()), points.get((j)%points.size())) +
            		  Evaluator.angle(points.get((i)%points.size()), points.get((j)%points.size()), points.get((j-1+points.size())%points.size())) +
            		  Evaluator.angle(points.get((i+1+points.size())%points.size()), points.get((j+1)%points.size()), points.get((j+2)%points.size())) +
            		  Evaluator.angle(points.get((j+1)%points.size()), points.get((i+1)%points.size()), points.get((i+2)%points.size()));
              
              
              if (a+b+a1>c+d+a2) {
            	  if(c < edgeThreshold && d < edgeThreshold) {
            		  ArrayList<Point> p=new ArrayList<Point>();
                      for (int k=0;k<=i;k++) p.add(points.get(k));
                      for (int k=j;k>i;k--) p.add(points.get(k));
                      for (int k=j+1;k<points.size();k++) p.add(points.get(k));
                      return p;
            	  }else {
            		  ArrayList<Point> p=new ArrayList<Point>();
            		  for (int k=0;k<=i;k++) p.add(points.get(k));
            		  if(points.get(i).distance(points.get(j)) > edgeThreshold) {
            			  ArrayList<Point> pts = expandPaths(original, points.get(i), points.get(j), paths, dist);
            			  pts.remove(0);
            			  p.addAll(pts);
            		  }
            		  
            		  for (int k=j;k>i;k--) p.add(points.get(k));
            		  
            		  if(points.get((i+1)%points.size()).distance(points.get((j+1)%points.size())) > edgeThreshold) {
            			  ArrayList<Point> pts = expandPaths(original, points.get((i+1)%points.size()), points.get((j+1)%points.size()), paths, dist);
            			  pts.remove(0);
            			  //Collections.reverse(pts);
            			  p.addAll(pts);
            		  }
            		  for (int k=j+1;k<points.size();k++) p.add(points.get(k));
                      return p;
            	  }
              }
          }
      }
      return points;
  }
  
  
  
  
  
  
  
  
  /*
  private double score(ArrayList<Point> points){
      return scoreDistance(points)+scoreAngle(points);
  }
  
  private double scoreDistance(ArrayList<Point> points) {
	  double s=0;
      for (int i=0;i<points.size();i++) s+=points.get(i).distance(points.get((i+1)%points.size()));
      return s;
  }
  
  private double scoreAngle(ArrayList<Point> points) {
	  double s=0;
	  double scal;
	  Point p,q,r;
	  FloatPoint v1, v2;
	  double value;
	  for(int i = 1; i <= points.size(); i++) {
		  p = points.get((i-1)%points.size());
		  q = points.get((i)%points.size());
		  r = points.get((i+1)%points.size());
		  v1 = new FloatPoint(q.getX()-p.getX(),q.getY()-p.getY());
		  v2 = new FloatPoint(r.getX()-q.getX(),r.getY()-q.getY());
		  scal = v1.x*v2.x + v1.y*v2.y;
		  value = scal/(v1.dist()*v2.dist());
		  if(value < -1) {
			  value = -1;
		  }else if(value > 1) {
			  value = 1;
		  }
		  s += (100./Math.PI)*Math.acos(value);
		  
	  }
	  return s;
  }*/
}

class FloatPoint{
	double x;
	double y;
	
	public FloatPoint(double x, double y) {
		this.x = x;
		this.y = y;
	}
	
	public double dist() {
		return Math.sqrt(x*x+y*y);
	}
	
}