//This scilab script builds the conductivity file for tutorial 3.2
//in the E4D users guide. 


//make sure there is sufficient memory in the stack
stacksize('max');


//___________________________READ THE MESH_____________________________________
//open the node file for writing
f1=mopen('mbsl.1.node','r');

//read the number of nodes from the first line of the node file
//and discard the rest of the line
nn=mfscanf(f1,'%f');              //nn is the number of nodes
mgetl(f1,1);                      //discard the rest, go to the second line
nods=zeros(nn,3);                 //allocate an array to hold node locations
bnds=zeros(nn,1);                 //allocate an array to hold node boundaries

//read the node locations and boundary values into the arrays
//columns 1 (j1) and 5 (j2) are not needed. 
//ns is a parameter that holds the number of elements successfully read 
for i=1:nn
 [ns,j1,nods(i,1),nods(i,2),nods(i,3),j2,bnds(i)]= ...
 mfscanf(f1,'%f %f %f %f %f %f');
end;
mclose(f1);

//Read the mesh translation values and translate the node locations
trn=read('mbsl.trn',1,3);
nods(:,1)=nods(:,1)+trn(1);
nods(:,2)=nods(:,2)+trn(2);
nods(:,3)=nods(:,3)+trn(3);

//open the element file
f1=mopen('mbsl.1.ele','r');

//read the number of elements from the first line
ne=mfscanf(f1,'%f');
//discard the rest of the first line
mgetl(f1,1);

//Allocate arrays for the elements and element zones.
//An element is defined by 4 nodes. Each row of the element
//file lists the 4 nodes defining that element.
el=zeros(ne,4);
zn=zeros(ne,1);

//Read the elements and zones
for i=1:ne
 [ns,j1,el(i,1),el(i,2),el(i,3),el(i,4),zn(i)]=...
 mfscanf(f1,'%f %f %f %f %f %f');
end;
mclose(f1);
//_____________________________________________________________________________


//_________________COMPUTE AND STORE THE ELEMENT MIDPOINTS_____________________
//Allocate an array to store the element centroids. Columns 1-3 are the
//x,y,and z positions of the element centroid.
mids=zeros(ne,3);

//Compute the centroid for each element. The centroid is the average
//position of the 4 nodes making the element. 
for i=1:ne
  mids(i,1)=.25*sum(nods(el(i,1:4),1));
  mids(i,2)=.25*sum(nods(el(i,1:4),2));
  mids(i,3)=.25*sum(nods(el(i,1:4),3));
end;
//_____________________________________________________________________________

//_______________________ASSIGN THE CONDUCTIVITIES_____________________________
//allocate an array to hold the conductivity values and assign the baseline
//value of 0.001 S/m
sig=0.001+zeros(ne,1);

//Check each element to determine if the centroid is within the 
//dipping plane, and if so set the conductivity to 0.01 S/m.
for i=1:ne

  //zt is the top of the dipping plane at this centroids x-value
  zt=.25*mids(i,1) - 1.75;

  //zb is the bottom of the dipping plane at this centroids x-value
  zb=.25*mids(i,1) - 2.25;

  //if the z-value of this centroid is between the top and bottom
  //of the dipping plane set the value to 0.01 S/m.
  if(mids(i,3)<=zt & mids(i,3)>=zb) then
     sig(i)=0.01;
  end;
end;
//_____________________________________________________________________________


//_________________________WRITE THE CONDUCTIVITY FILE_________________________
f1=mopen('mbsl_truesig.sig','w');
mfprintf(f1,'%-10.0f 1\n',ne);
for i=1:ne
 mfprintf(f1,'%-5.4f\n',sig(i));
end;
mclose(f1);
//_____________________________________________________________________________

