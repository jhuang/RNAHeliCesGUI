package fr.orsay.lri.varna.models.templates;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import fr.orsay.lri.varna.exceptions.ExceptionInvalidRNATemplate;
import fr.orsay.lri.varna.models.rna.RNA;
import fr.orsay.lri.varna.models.templates.RNATemplate.RNATemplateElement;
import fr.orsay.lri.varna.models.templates.RNATemplate.RNATemplateHelix;
import fr.orsay.lri.varna.models.templates.RNATemplate.RNATemplateUnpairedSequence;
import fr.orsay.lri.varna.models.treealign.AlignedNode;
import fr.orsay.lri.varna.models.treealign.RNANodeValue;
import fr.orsay.lri.varna.models.treealign.RNANodeValue2;
import fr.orsay.lri.varna.models.treealign.RNATree2;
import fr.orsay.lri.varna.models.treealign.RNATree2Exception;
import fr.orsay.lri.varna.models.treealign.Tree;
import fr.orsay.lri.varna.models.treealign.TreeAlign;
import fr.orsay.lri.varna.models.treealign.TreeAlignException;
import fr.orsay.lri.varna.models.treealign.TreeAlignResult;

/**
 * This class is about the alignment between a tree of RNANodeValue2
 * and a tree of RNANodeValueTemplate.
 * 
 * @author Raphael Champeimont
 */
public class RNATemplateAlign {
	
	private static boolean canBePartOfAnHelix(RNANodeValue2 leftNodeValue) {
		return (leftNodeValue != null) && leftNodeValue.isSingleNode();
	}
	
	private static boolean canBePartOfASequence(RNANodeValue2 leftNodeValue) {
		return (leftNodeValue != null) && !leftNodeValue.isSingleNode();
	}
	
	
	/**
	 * This method takes an alignment between a tree of RNANodeValue2
	 * of RNANodeValue and a tree of RNANodeValue2 of RNANodeValueTemplate.
	 * It returns the corresponding RNATemplateMapping.
	 */
	public static RNATemplateMapping makeTemplateMapping(Tree<AlignedNode<RNANodeValue2,RNANodeValueTemplate>> alignment) throws RNATemplateMappingException {
		RNATemplateMapping mapping = new RNATemplateMapping();
		
		RNATemplateHelix currentHelix = null;
		LinkedList<Tree<AlignedNode<RNANodeValue2,RNANodeValueTemplate>>> remainingNodes = new LinkedList<Tree<AlignedNode<RNANodeValue2,RNANodeValueTemplate>>>();
		// The reason why this algorithm is not trivial is that we may have
		// a longer helix on the RNA side than on the template side, in which
		// case some nodes on the RNA side are going to be alone while we
		// would want them to be part of the helix.
		List<RNANodeValue2> nodesInSameHelix = new LinkedList<RNANodeValue2>();
		remainingNodes.add(alignment);
		while (!remainingNodes.isEmpty()) {
			Tree<AlignedNode<RNANodeValue2,RNANodeValueTemplate>> node = remainingNodes.getLast();
			remainingNodes.removeLast();
			
			Tree<RNANodeValue2> leftNode = node.getValue().getLeftNode();
			Tree<RNANodeValueTemplate> rightNode = node.getValue().getRightNode();
			if (leftNode != null && leftNode.getValue() != null) {
				RNANodeValue2 leftNodeValue = leftNode.getValue();
				// We have a real left (RNA side) node
				
				if (rightNode != null && rightNode.getValue() != null) {
					// We have a real right (template side) node
					RNANodeValueTemplate rightNodeValue = rightNode.getValue();
					
					if (rightNodeValue instanceof RNANodeValueTemplateBasePair
							&& canBePartOfAnHelix(leftNodeValue)) {
						RNATemplateHelix helix = ((RNANodeValueTemplateBasePair) rightNodeValue).getHelix();
						currentHelix = helix;
						mapping.addCouple(leftNodeValue.getNode().getLeftBasePosition(), helix);
						mapping.addCouple(leftNodeValue.getNode().getRightBasePosition(), helix);
						
						// Maybe we have marked nodes as part of the same helix
						// when we didn't know yet which helix it was.
						for (RNANodeValue2 v: nodesInSameHelix) {
							mapping.addCouple(v.getNode().getLeftBasePosition(), helix);
							mapping.addCouple(v.getNode().getRightBasePosition(), helix);
						}
						nodesInSameHelix.clear();
					} else if (rightNodeValue instanceof RNANodeValueTemplateSequence
							&& canBePartOfASequence(leftNodeValue)) {
						currentHelix = null;
						nodesInSameHelix.clear();
						RNATemplateUnpairedSequence sequence = ((RNANodeValueTemplateSequence) rightNodeValue).getSequence();
						for (RNANodeValue nodeValue: leftNode.getValue().getNodes()) {
							mapping.addCouple(nodeValue.getLeftBasePosition(), sequence);
						}
					} else {
						currentHelix = null;
						nodesInSameHelix.clear();
					}
				} else {
					// We have no right (template side) node
					
					if (canBePartOfAnHelix(leftNodeValue)) {
						if (currentHelix != null) {
							// We may be in this case if the RNA sequence
							// contains a longer helix than in the template
							mapping.addCouple(leftNodeValue.getNode().getLeftBasePosition(), currentHelix);
							mapping.addCouple(leftNodeValue.getNode().getRightBasePosition(), currentHelix);
						} else {
							// Maybe this left node is part of an helix
							// which is smaller in the template
							nodesInSameHelix.add(leftNodeValue);
						}
					} else {
						currentHelix = null;
						nodesInSameHelix.clear();
					}
				}
			} else {
				currentHelix = null;
				nodesInSameHelix.clear();
			}
			
			// If this node has children, add them in the stack
			List<Tree<AlignedNode<RNANodeValue2,RNANodeValueTemplate>>> children = node.getChildren();
			int n = children.size();
			if (n > 0) {
				for (int i=n-1; i>=0; i--) {
					// We add the children in their reverse order so they
					// are given in the original order by the iterator
					remainingNodes.add(children.get(i));
				}
			} else {
				// We will now stop going down, so there is no reason
				// the current helix will still be the same
				currentHelix = null;
				nodesInSameHelix.clear();
			}
			
		}
		
		return mapping;
	}
	
	
	
	public static void printMapping(RNATemplateMapping mapping, RNATemplate template, String sequence) {
		Iterator<RNATemplateElement> iter = template.rnaIterator();
		while (iter.hasNext()) {
			RNATemplateElement element = iter.next();
			System.out.println(element.toString());
			ArrayList<Integer> A = mapping.getAncestor(element);
			if (A != null) {
				RNATemplateAlign.printIntArrayList(A);
				for (int n=A.size(), i=0; i<n; i++) {
					System.out.print("\t" + sequence.charAt(A.get(i)));
				}
				System.out.println("");
			} else {
				System.out.println("\tno match");
			}
		}
	}
	
	
	/**
	 * Align the given RNA with the given RNA template.
	 */
	public static TreeAlignResult<RNANodeValue2,RNANodeValueTemplate> alignRNAWithTemplate(RNA rna, RNATemplate template) throws RNATemplateDrawingAlgorithmException {
		try {
			Tree<RNANodeValue2> rnaAsTree = RNATree2.RNATree2FromRNA(rna);
			Tree<RNANodeValueTemplate> templateAsTree = template.toTree();
			TreeAlign<RNANodeValue2,RNANodeValueTemplate> treeAlign = new TreeAlign<RNANodeValue2,RNANodeValueTemplate>(new RNANodeValue2TemplateDistance());
			TreeAlignResult<RNANodeValue2,RNANodeValueTemplate> result = treeAlign.align(rnaAsTree, templateAsTree);
			return result;
		} catch (RNATree2Exception e) {
			throw (new RNATemplateDrawingAlgorithmException("RNATree2Exception: " + e.getMessage()));
		} catch (ExceptionInvalidRNATemplate e) {
			throw (new RNATemplateDrawingAlgorithmException("ExceptionInvalidRNATemplate: " + e.getMessage()));
		} catch (TreeAlignException e) {
			throw (new RNATemplateDrawingAlgorithmException("TreeAlignException: " + e.getMessage()));
		}
	}
	
	/**
	 * Build a mapping between an RNA sequence and a template using
	 * a single tree alignment.
	 */
	public static RNATemplateMapping singlePassMap(RNA rna, RNATemplate template) throws RNATemplateDrawingAlgorithmException {
		try {
			TreeAlignResult<RNANodeValue2,RNANodeValueTemplate> alignResult = RNATemplateAlign.alignRNAWithTemplate(rna, template);
			Tree<AlignedNode<RNANodeValue2,RNANodeValueTemplate>> alignment = alignResult.getAlignment();
			RNATemplateMapping mapping = RNATemplateAlign.makeTemplateMapping(alignment);
			return mapping;
		} catch (RNATemplateMappingException e) {
			throw (new RNATemplateDrawingAlgorithmException("RNATemplateMappingException: " + e.getMessage()));
		}
	}
	
	/**
	 * Build a mapping between an RNA sequence and a template using several
	 * tree alignments in order to support pseudoknots.
	 */
	public static RNATemplateMapping multiPassMap(RNA rna, RNATemplate template) throws RNATemplateDrawingAlgorithmException {
		return singlePassMap(rna, template);
		// TODO
	}
	
	
	

	/**
	 * Print an integer array.
	 */
	public static void printIntArray(int[] A) {
		for (int i=0; i<A.length; i++) {
			System.out.print("\t" + A[i]);
		}
		System.out.println("");
	}

	/**
	 * Print an integer ArrayList.
	 */
	public static void printIntArrayList(ArrayList<Integer> A) {
		for (int i=0; i<A.size(); i++) {
			System.out.print("\t" + A.get(i));
		}
		System.out.println("");
	}

	/**
	 * Print an matrix of shorts.
	 */
	public static void printShortMatrix(short[][] M) {
		System.out.println("Begin matrix");
		for (int i=0; i<M.length; i++) {
			for (int j=0; j<M[i].length; j++) {
				System.out.print("\t" + M[i][j]);
			}
			System.out.println("");
		}
		System.out.println("End matrix");
	}
	
	/**
	 * Convert a list of integers into an array of integers.
	 * The returned arrays is freshly allocated.
	 * Returns null if given null.
	 */
	public static int[] intArrayFromList(List<Integer> l) {
		if (l != null) {
			int n = l.size();
			int[] result = new int[n];
			for (int i=0; i<n; i++) {
				result[i] = l.get(i);
			}
			return result;
		} else {
			return null;
		}
	}
	
}
