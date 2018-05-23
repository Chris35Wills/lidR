#include "BoundingBox3D.h"
#include "Point.h"
#include "QuadTree3D.h"
#include "TreeCollection.h"
#include "Tree.h"

using namespace Rcpp;

typedef boost::geometry::model::point<double, 2, boost::geometry::cs::cartesian> point_t;

//========================================================================================
template<typename T> void apply2DFilter( std::vector<T> &subProfile, std::vector<T> &subProfileSubset );
TreeCollection<PointXYZ> PTrees_segmentation( S4 &las, int k );
IntegerMatrix createCombination(int N);
template<typename T> double getScoreCombination( Tree<T> &t1, Tree<T> &t2, int k );
template<typename T> Tree<T> mergeTree( Tree<T> &t1, Tree<T> &t2, int k );
//========================================================================================


struct Increment {
  Increment() : m_value( 1 ) { }
  int operator()() { return m_value++; }
  int m_value;
};
//========================================================================================

// [[Rcpp::export]]
std::vector<int> lastrees_PTrees( S4 las, NumericVector k_values  )
{
  std::vector<std::vector<double> > M;
  //sorting k values in decreasing order
  std::sort(k_values.begin(), k_values.end(), std::greater<double>());
  Rcpp::Rcout << "k = "<< k_values[0] << std::endl;
  TreeCollection<PointXYZ> trees_kRef = PTrees_segmentation( las, k_values[0]);

  //Rcpp::Rcout << "===== END SEGMENTATION 1 =====" << std::endl;

  double total = 0;
  for (int i = 0; i < trees_kRef.nbTree; i++)
    total += trees_kRef.treeStorage[i].nbPoints;

  std::vector<int> idResult( total,0);
  for (int nb_k = 1; nb_k < k_values.size(); nb_k++)
  {
    int index = 1;
    Rcpp::Rcout << "k = "<< k_values[nb_k] << std::endl;
    TreeCollection<PointXYZ> trees_kToCompare = PTrees_segmentation( las, k_values[nb_k]);

    Rcpp::Rcout << "===== END SEGMENTATION =====" << std::endl;
    std::vector<double> treeScores;
    std::vector<int> treeIDs;
    std::vector<PointXYZ> treeZmax;
    TreeCollection<PointXYZ> trees_kResult;
    //for (int t = 0; t < trees_kRef.nbTree; t++)
    Rcpp::Rcout << "trees_kRef.nbTree="<<trees_kRef.nbTree<<std::endl;
    for (int t = 0; t < trees_kRef.nbTree; t++)
    {
      M.clear();
      Rcpp::Rcout <<"=========="<<std::endl;
      treeScores.clear();
      treeIDs.clear();
      treeZmax.clear();
      polygon poly2D;

      //Création du polygone de recherche poly2D de l'arbre de référence
      //Cas où l'arbre n'est constitué que d'un ou deux points et n'a donc pas d'enveloppe convexe
      //On ajoute ce(s) point(s) comme un arbre isolé (-->?)
      if (trees_kRef.treeStorage[t].pointsCH.size() < 3)
      {
        //Rcpp::Rcout << "No pointsCH" << std::endl;
        trees_kResult.addTree( trees_kRef.treeStorage[t] );
        //Attribuer à tous les ID des points de cet arbre une meme valeur
        trees_kRef.treeStorage[t].editIdResult (idResult, index);
      }
      else //Cas où l'arbre de référence possède une enveloppe convexe
      {
        boost::geometry::assign_points(poly2D, trees_kRef.treeStorage[t].pointsCH);
        point_t p;

        //Dans la boucle qui suit, on regarde ensuite quels apex de la collection d'arbres à comparer
        //se situent dans le polygone défini précédemment
        int nbApex = 0;
        for (int n = 0; n < trees_kToCompare.nbTree; n++)
        {
          //Apex search for each tree
          PointXYZ apex;
          trees_kToCompare.treeStorage[n].getZMax( apex );

          //Transformation des apex de pointXYZ en boost::points
          boost::geometry::set<0>(p, apex.x);
          boost::geometry::set<1>(p, apex.y);
          //Apex inclus dans polygone?
          std::vector<double> tmp;
          if (boost::geometry::covered_by(p, poly2D)  == TRUE)
          {
            treeScores.push_back(trees_kToCompare.treeStorage[n].scoreGlobal);  //stockage du score global de l'arbre
            treeIDs.push_back(n);               //stockage de l'ID correspondant
            treeZmax.push_back(trees_kToCompare.treeStorage[n].Zmax);           //stockage de la hauteur max
            nbApex ++;
            tmp.clear();
            tmp.push_back(boost::geometry::get<0>(p));
            tmp.push_back(boost::geometry::get<1>(p));
            M.push_back(tmp);
          }
        }
        //Rcpp::Rcout << "numberTree=" << nbApex<< std::endl;
        //Rcpp::Rcout << "treeIDs = ";
        /*for (int i = 0; i < treeIDs.size(); i++)
          Rcpp::Rcout << treeIDs[i]<<" /";
        Rcpp::Rcout << std::endl;*/


        double scoreRef = trees_kRef.treeStorage[t].scoreGlobal;
        //Rcpp::Rcout << "scoreRef= " << scoreRef << std::endl;
        //If one apex of tree n is found inside convex hull of tree t (page 103 Fig5 PartA-1)
        if(treeIDs.size() == 1)
        {
          //Rcpp::Rcout << "Only one apex" << std::endl;
          int keep_id = treeIDs[0];
          double scoreToCompare = trees_kToCompare.treeStorage[keep_id].scoreGlobal;
          //Rcpp::Rcout << "scoreToCompare= " << scoreToCompare<<" "<<treeScores[0]<<std::endl;
          //Comparison of score A and score D
          if (scoreRef > scoreToCompare)
          {
            //Rcpp::Rcout << "Cas1 " << std::endl;
            trees_kResult.addTree( trees_kRef.treeStorage[t] );
            //Attribuer à tous les ID des points de cet arbre une meme valeur
            trees_kRef.treeStorage[t].editIdResult (idResult, index);
          }
          else
          {
            //Rcpp::Rcout << "Cas2 " << std::endl;
            trees_kResult.addTree( trees_kToCompare.treeStorage[keep_id] );
            //Attribuer à tous les ID des points de cet arbre une meme valeur
            trees_kToCompare.treeStorage[keep_id].editIdResult (idResult, index);
          }
        }
        //If two apices of tree n were found inside convex hull of tree t (page 103 Fig5 PartA-2)
        else if (treeIDs.size() == 2)
        {
          //Rcpp::Rcout << "Only two apex" << std::endl;
          //Comparison of score B and score E+J
          if ( (treeScores[0] + treeScores[1])/2.0 > scoreRef)
          {
            //Rcpp::Rcout << "treeScores[0] + treeScores[1])/2 =" << (treeScores[0] + treeScores[1])/2.0 << std::endl;
            trees_kResult.addTree( trees_kToCompare.treeStorage[treeIDs[0]] );
            trees_kResult.addTree( trees_kToCompare.treeStorage[treeIDs[1]] );
            //Attribuer à tous les ID des points de cet arbre une meme valeur
            trees_kToCompare.treeStorage[treeIDs[0]].editIdResult (idResult, index);
            trees_kToCompare.treeStorage[treeIDs[1]].editIdResult (idResult, index);
          }
          else
          {
            trees_kResult.addTree( trees_kRef.treeStorage[t] );
            //Attribuer à tous les ID des points de cet arbre une meme valeur
            trees_kRef.treeStorage[t].editIdResult (idResult, index);
          }
        }
        //If more than two apices of tree n were found inside convex hull of tree t
        //--> test of each combination
        else if (treeIDs.size() > 2)
        {
          //Rcpp::Rcout << "More than two apex" << std::endl;
          Rcpp::IntegerMatrix combination = createCombination( treeIDs.size() );

          //Moyenne totale des scores de l'arbre à comparer
          double sumAllTrees = std::accumulate(treeScores.begin(), treeScores.end(), 0.0);
          double maxScore = sumAllTrees / treeScores.size();
          //Rcpp::Rcout << "maxScore=" <<maxScore <<std::endl;

          std::vector<int> includedTreeForScore (2, INT16_MIN);
          std::vector<int> remainingTreeIDs(treeIDs.size(), 0);
          //Rcpp::Rcout << "===1===" << std::endl;

          //Recherche de toutes les combinaisons d'arbres possibles et calcul de leurs scores
          //le meilleur score est conservé dans "maxScore" et sa combinaison associée est stockée dans "includedTreeForScore"
          for (int i = 0; i < combination.nrow(); i++)
          {
            //Rcpp::Rcout << i <<" ";
            if (combination(i,0) != 0)
            {
              int nbTreeConsidered = 0;
              double calculatedScoreForCombination = 0;
              int tree1 = treeIDs[combination(i,0)-1];
              int tree2 = treeIDs[combination(i,1)-1];
              //Rcpp::Rcout << "tree1 =" << tree1 << "tree2="<<tree2 <<std::endl;
              double scoreTreeCombination = getScoreCombination( trees_kToCompare.treeStorage[tree1], trees_kToCompare.treeStorage[tree2], nb_k );

              //Rcpp::Rcout << scoreTreeCombination <<" "<<std::endl;
              //ajouter à scoreTreeCombination la moyenne des restants --> finalCombinationScore
              remainingTreeIDs.clear();
              remainingTreeIDs.assign(treeIDs.size(), 0);
              std::generate(remainingTreeIDs.begin(), remainingTreeIDs.end(), Increment());
              remainingTreeIDs.erase(std::remove(remainingTreeIDs.begin(), remainingTreeIDs.end(), combination(i,0)), remainingTreeIDs.end());
              remainingTreeIDs.erase(std::remove(remainingTreeIDs.begin(), remainingTreeIDs.end(), combination(i,1)), remainingTreeIDs.end());

              double otherTreeScore = 0;
              for (int j = 0; j < remainingTreeIDs.size(); j++)
                otherTreeScore += treeScores[remainingTreeIDs[j]-1];

              double finalCombinationScore = ((otherTreeScore / remainingTreeIDs.size()) + scoreTreeCombination) / 2.0;
              //Rcpp::Rcout << "finalCombinationScore=" << finalCombinationScore  << std::endl;

              if (maxScore < finalCombinationScore )
              {
                maxScore = finalCombinationScore;
                for (int k = 0; k < includedTreeForScore.size(); k++)
                  includedTreeForScore[k] = combination(i,k);
              }
            }
          }

          //Comparison between max score of possible combination (including mean score of all individual trees) and reference score
          if (maxScore <= scoreRef) //storage of tree coming from reference tree collection (the one which creates polygon for apex search)
          {
            trees_kResult.addTree( trees_kRef.treeStorage[t] );
            trees_kRef.treeStorage[t].editIdResult (idResult, index);
          }
          else  //storage of best combination or all tree coming from tree collection to compare
          {
            //Si le maxScore est le même que celui de la moyenne de tous les scores des arbres à comparer
            //et qu'aucune combinaison n'a été stockée --> on stocke tous les arbres isolés
            if (maxScore == (sumAllTrees / treeScores.size()) && includedTreeForScore[0] == INT16_MIN )
            {
              for (int i = 0; i < treeIDs.size(); i++)
              {
                trees_kResult.addTree( trees_kToCompare.treeStorage[treeIDs[i]] );
                trees_kToCompare.treeStorage[treeIDs[i]].editIdResult(idResult, index);
              }
            }
            else  //sinon on stocke la combinaison (deux arbres fusionnés) et tous les autres arbres isolés
            {
              //Ajout de l'arbre fusionné (une des combinaisons de deux arbres isolés)
              int tree1 = treeIDs[includedTreeForScore[0]-1];
              int tree2 = treeIDs[includedTreeForScore[1]-1];
              Tree<PointXYZ> combinedTree;
              combinedTree = mergeTree(trees_kToCompare.treeStorage[tree1], trees_kToCompare.treeStorage[tree2], nb_k);

              trees_kResult.addTree( combinedTree );
              combinedTree.editIdResult (idResult, index);

              //ajout des arbres restants isolés
              std::vector<int> remainingTreeIDs(treeIDs.size());
              std::generate(remainingTreeIDs.begin(), remainingTreeIDs.end(), Increment());
              remainingTreeIDs.erase(std::remove(remainingTreeIDs.begin(), remainingTreeIDs.end(), includedTreeForScore[0]), remainingTreeIDs.end());
              remainingTreeIDs.erase(std::remove(remainingTreeIDs.begin(), remainingTreeIDs.end(), includedTreeForScore[1]), remainingTreeIDs.end());
              for (int i = 0; i < remainingTreeIDs.size(); i++)
              {
                int ID = remainingTreeIDs[i] - 1;
                trees_kResult.addTree( trees_kToCompare.treeStorage[treeIDs[ID]] );
                trees_kToCompare.treeStorage[treeIDs[ID]].editIdResult (idResult, index);
              }
            }
          }
        }
      }

      for (int i = 0; i < trees_kResult.idTreeStorage.size(); i++)
      {
        idResult.push_back(trees_kResult.idTreeStorage[i]);
      }
    }
  }
  return(idResult);
}

//========================================================================================
//                               PTREES_SEGMENTATION
//========================================================================================
TreeCollection<PointXYZ> PTrees_segmentation( S4 &las, int k )
{
  //Data conversion from las object to vector of PointXYZ
  DataFrame data = as<Rcpp::DataFrame>(las.slot("data"));
  NumericVector X = data["X"];
  NumericVector Y = data["Y"];
  NumericVector Z = data["Z"];

  std::vector<PointXYZ> points(X.size());
  for (int i = 0 ; i < X.size() ; i++)
    points[i] = PointXYZ(X[i], Y[i], Z[i], i);

  //Vector sorting by Z (from max to min)
  std::sort( points.begin(), points.end(), ZSortPointBis<PointXYZ>() );

  //First Zmax defines first tree in trees
  PointXYZ pointToSort = points[0];
  Tree<PointXYZ> treeInit( pointToSort );
  TreeCollection<PointXYZ> trees ( treeInit );

  //Creation of vector that stores relation between tree number and TreeCollection
  //Update for first point assignation
  std::vector<int> idTree( X.size(), 0 );
  idTree[pointToSort.id] = trees.nbTree;
  Rcpp::Rcout << "NumeroArbre="<< trees.nbTree <<" ";

  //Creation of a QuadTree
  QuadTree3D<PointXYZ> *treeOI;
  treeOI = QuadTreeCreate(points);

  //Setting of resolution for QuadTree--> put as argument? resolution to use?
  PointXYZ reso = PointXYZ(0.01,0.01,0.01);

  std::vector<int> knnTreeID;
  std::vector<PointXYZ> result, filteredResult;
  BoundingBox3D<PointXYZ> bbox;

  //Each point is attributed to an existing or new tree
  int limit = points.size();
  for (int i = 1; i < limit; i ++ )
  {
    //Initialisation of point to assign in trees
    pointToSort = points[i];

    //Searching for k-nearest neighbours and storage of corresponding points into result
    bbox = BoundingBox3D<PointXYZ> (pointToSort, reso);
    result.clear();
    filteredResult.clear();
    treeOI->knn_lookup3D( pointToSort, k, result );  //result contains founded neighbours + point of interest

    //Removal of points having a planimetric distance from pointToSort above threshold T (page 100 Eq. 1/2)
    apply2DFilter(result, filteredResult);

    //Searching if these k points were already classified in idTree using ID number in result
    knnTreeID.clear();
    knnTreeID.assign( filteredResult.size(), 0 );
    for (int n = 0; n < filteredResult.size(); n++ )
    {
      knnTreeID[n] = idTree[filteredResult[n].id];
    }

    //Removal of duplicates tree IDs and index value 0
    sort( knnTreeID.begin(), knnTreeID.end() );
    knnTreeID.erase( unique( knnTreeID.begin(), knnTreeID.end() ), knnTreeID.end() );
    knnTreeID.erase(knnTreeID.begin());

    //Three possibilities for next classification:
    //1) If no classified points are found in the neighbourhood (page 101 Fig4.B situation 1)
    if ( knnTreeID.empty() == TRUE )
    {
      Tree<PointXYZ> newTree( pointToSort );
      trees.addTree( newTree );
      idTree[pointToSort.id] = trees.nbTree;
    }
    //2) If only one identified tree in the neighbourhood (page 101 Fig4.B situation 2)
    else if ( knnTreeID.size() == 1)
    {
      trees.updateTree( knnTreeID[0], pointToSort );
      idTree[pointToSort.id] = knnTreeID[0];
    }
    //3) If several identified trees in the neighbourhood (page 101 Fig4.B situation 3)
    else
    {
      //To know which search Method is required for this tree subset (distance or area evaluation)
      //scan of tree numbers to identify if there is at least one with less than 2 points
      int searchMethod = 1;
      for ( int i = 0; i < knnTreeID.size(); i++ )
      {
        if ( trees.treeStorage[knnTreeID[i]-1].nbPoints < 2 )
          searchMethod = 2;
      }

      //Depending on previous result, selection of adapted searchMethod
      int resultID = knnTreeID[0];
      double areaValue = 0, distValue = 0;
      double diffHeight = 0;
      double thresholdZ = 5;                  //page 100 last paragraph
      switch(searchMethod)
      {
      case 1: { boost::geometry::model::ring<point_t> hull;
        trees.searchID_usingArea( knnTreeID, pointToSort, resultID, areaValue, hull );

        //Before association of pointToSort to best tree result,
        //testing if Z difference between pointToSort and lowest point in tree is under
        //a height difference threshold fixed at 5m --> page 100 last paragraph
        diffHeight = fabs(trees.treeStorage[resultID-1].findZMin() - pointToSort.z);
        if ( diffHeight <= thresholdZ )
        {
          trees.treeStorage[resultID-1].addPoint( pointToSort, areaValue, hull );
          trees.individualTreeSize[resultID-1]++;
          idTree[pointToSort.id] = resultID;
        }
        else
        {
          Tree<PointXYZ> newTree( pointToSort );
          trees.addTree( newTree );
          idTree[pointToSort.id] = trees.nbTree;
        }
        break;
      }

      case 2: { trees.searchID_usingDist( knnTreeID, pointToSort, resultID, distValue );

        //Before association of pointToSort to best tree result,
        //testing if Z difference between pointToSort and lowest point in tree is under
        //a height difference threshold fixed at 5m --> page 100 last paragraph
        diffHeight = std::abs(trees.treeStorage[resultID-1].findZMin() - pointToSort.z);
        if ( diffHeight <= thresholdZ )
        {
          trees.treeStorage[resultID-1].addPoint_dist( pointToSort, distValue );
          idTree[pointToSort.id] = resultID;
          trees.individualTreeSize[resultID-1]++;
        }
        else
        {
          Tree<PointXYZ> newTree( pointToSort );
          trees.addTree( newTree );
          idTree[pointToSort.id] = trees.nbTree;
        }
        break;
      }
      }

    }

  }
  trees.idTreeStorage = idTree;
  trees.calculateTreeScores(k);

  delete treeOI;

return(trees);
}



//========================================================================================
//                              APPLY 2D FILTER FUNCTION
//========================================================================================
//Function that calculates 2D distance between one point (first one in subProfile) and its nearest neighbours
//returns only thoose under a calculated threshold
template<typename T> void apply2DFilter( std::vector<T> &subProfile, std::vector<T> &subProfileSubset )
{
  double meanValueForThreshold = 0;
  std::vector<double> dist;
  double stdValueForThreshold = 0;
  double val = 0;
  //Euclidian Distance calculation in 2D for all neighbours regarding the reference point (storage in Z value)
  for ( int i = 1; i < subProfile.size(); i++ )
  {
    val = euclidianDistance2D_inZ( subProfile[0], subProfile[i] );
    meanValueForThreshold += val;
    dist.push_back( val );
  }

  //--------------------------------------------------------------------------------------
  //Threshold definition (page 100 'segmentation principles')
  //defined as the mean plus twice the std of the planimetric distances of a subset of points
  meanValueForThreshold /= (double)(subProfile.size() - 1);
  double sum = 0;
  for(int i = 0; i<dist.size(); i++)
    sum += (dist[i]-meanValueForThreshold) * (dist[i]-meanValueForThreshold);

  stdValueForThreshold = sqrt( sum / (double) dist.size() );

  double threshold = meanValueForThreshold + 2*stdValueForThreshold;
  //--------------------------------------------------------------------------------------
  //Keeping all points below threshold and storage of their IDs
  int i = 0, keep = 0;
  subProfileSubset.push_back(subProfile[0]);
  while ( (i < subProfile.size()-1) && (dist[i] <= threshold) )
  {
    subProfileSubset.push_back( subProfile[i+1] );
    i++;
  }

}

template<typename T> Tree<T> mergeTree( Tree<T> &t1, Tree<T> &t2, int k )
{
  Tree<T> newTree(t1);
  for (int i = 0; i < t2.nbPoints; i++)
    newTree.points.push_back(t2.points[i]);

  newTree.updateArea();
  newTree.getScore(k);
  newTree.nbPoints += t2.nbPoints;
  return(newTree);
}

//========================================================================================
template<typename T> double getScoreCombination( Tree<T> &t1, Tree<T> &t2, int k )
{
  Tree<T> combinedTrees;
  combinedTrees = mergeTree( t1, t2, k );
  return (combinedTrees.scoreGlobal);
}

//========================================================================================

//========================================================================================
//                              CREATE COMBINATION
//========================================================================================
//https://stackoverflow.com/questions/12991758/creating-all-possible-k-combinations-of-n-items-in-c
//derived from Rosetta Code
IntegerMatrix createCombination(int N)
{
  IntegerMatrix res(pow(N-1,2)-(N-2), 2);   //pas optimal
  int ind = 0, ind2 = 0;
  //for (int j = 1; j< N; j++)
  //{
    //ind = 0;
    std::string bitmask(2, 1); // K leading 1's
    bitmask.resize(N, 0); // N-K trailing 0's
    // print integers and permute bitmask
    do {
      for (int i = 0; i < N+1; i++) // [0..N-1] integers
      {
        if (bitmask[i]) res(ind,ind2++) = i+1;
      }
      ind++;
      ind2 = 0;
    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));
  //}
  return(res);
}




