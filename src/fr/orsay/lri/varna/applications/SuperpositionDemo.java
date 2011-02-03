/*
 VARNA is a tool for the automated drawing, visualization and annotation of the secondary structure of RNA, designed as a companion software for web servers and databases.
 Copyright (C) 2008  Kevin Darty, Alain Denise and Yann Ponty.
 electronic mail : Yann.Ponty@lri.fr
 paper mail : LRI, bat 490 Université Paris-Sud 91405 Orsay Cedex France

 This file is part of VARNA version 3.1.
 VARNA version 3.1 is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

 VARNA version 3.1 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along with VARNA version 3.1.
 If not, see http://www.gnu.org/licenses.
 */
package fr.orsay.lri.varna.applications;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

import fr.orsay.lri.varna.VARNAPanel;
import fr.orsay.lri.varna.exceptions.ExceptionDrawingAlgorithm;
import fr.orsay.lri.varna.exceptions.ExceptionFileFormatOrSyntax;
import fr.orsay.lri.varna.exceptions.ExceptionNAViewAlgorithm;
import fr.orsay.lri.varna.exceptions.ExceptionNonEqualLength;
import fr.orsay.lri.varna.exceptions.ExceptionUnmatchedClosingParentheses;
import fr.orsay.lri.varna.exceptions.MappingException;
import fr.orsay.lri.varna.models.rna.Mapping;
import fr.orsay.lri.varna.models.rna.RNA;

public class SuperpositionDemo extends JFrame {

	/**
	 * 
	 */
	private static final long serialVersionUID = -790155708306987257L;

	private static final String DEFAULT_SEQUENCE = "CAGCACGACACUAGCAGUCAGUGUCAGACUGCAIACAGCACGACACUAGCAGUCAGUGUCAGACUGCAIACAGCACGACACUAGCAGUCAGUGUCAGACUGCAIA";

	private static final String DEFAULT_STRUCTURE1 = "..(((((...(((((...(((((...(((((.....)))))...))))).....(((((...(((((.....)))))...))))).....)))))...)))))..";
	private static final String DEFAULT_STRUCTURE2 = "..(((((...(((((...(((((........(((((...(((((.....)))))...)))))..................))))).....)))))...)))))..";
	// private static final String DEFAULT_STRUCTURE1 = "((((....))))";
	// private static final String DEFAULT_STRUCTURE2 =
	// "((((..(((....)))..))))";

	private VARNAPanel _vp;

	private JPanel _tools = new JPanel();
	private JPanel _input = new JPanel();

	private JPanel _seqPanel = new JPanel();
	private JPanel _struct1Panel = new JPanel();
	private JPanel _struct2Panel = new JPanel();
	private JLabel _info = new JLabel();
	private JTextField _struct1 = new JTextField(DEFAULT_STRUCTURE1);
	private JTextField _struct2 = new JTextField(DEFAULT_STRUCTURE2);
	private JTextField _seq = new JTextField(DEFAULT_SEQUENCE);
	private JLabel _struct1Label = new JLabel(" Str1:");
	private JLabel _struct2Label = new JLabel(" Str2:");
	private JLabel _seqLabel = new JLabel(" Seq:");
	private JButton _goButton = new JButton("Go");
	private JButton _switchButton = new JButton("Switch");

	private String _str1Backup = "";
	private String _str2Backup = "";
	private RNA _RNA1 = new RNA();
	private RNA _RNA2 = new RNA();

	private static String errorOpt = "error";
	@SuppressWarnings("unused")
	private boolean _error;

	private Color _backgroundColor = Color.white;

	@SuppressWarnings("unused")
	private int _algoCode;

	private int _currentDisplay = 1;

	public SuperpositionDemo() {
		super();
		try {
			_vp = new VARNAPanel(_seq.getText(), getStruct1());
		} catch (ExceptionNonEqualLength e) {
			_vp.errorDialog(e);
		}
		_vp.setPreferredSize(new Dimension(400, 400));
		RNAPanelDemoInit();
	}

	private void RNAPanelDemoInit() {
		int marginTools = 40;

		setBackground(_backgroundColor);
		_vp.setBackground(_backgroundColor);

		try {
			_RNA1.setRNA(_seq.getText(), getStruct1());
			_RNA2.setRNA(_seq.getText(), getStruct2());

			_vp.drawRNA(getRNA1());
		} catch (ExceptionUnmatchedClosingParentheses e1) {
			_vp.errorDialog(e1);
		} catch (ExceptionFileFormatOrSyntax e1) {
			_vp.errorDialog(e1);
		}

		Font textFieldsFont = Font.decode("MonoSpaced-PLAIN-12");

		_seqLabel.setHorizontalTextPosition(JLabel.LEFT);
		_seqLabel.setPreferredSize(new Dimension(marginTools, 15));
		_seq.setFont(textFieldsFont);
		_seq.setText(getRNA1().getSeq());

		_goButton.addActionListener(new ActionListener() {

			public void actionPerformed(ActionEvent e) {
				_currentDisplay = (_currentDisplay + 1) % 2;
				if (_currentDisplay == 1) {
					_vp.drawRNA(getRNA1());
				} else {
					_vp.drawRNA(getRNA2());
				}
				_vp.repaint();
			}
		});

		_switchButton.addActionListener(new ActionListener() {

			public void actionPerformed(ActionEvent e) {
				try {
					_currentDisplay = (_currentDisplay + 1) % 2;

					if (_currentDisplay == 1) {
						try {
							Mapping m = Mapping.readMappingFromAlignment(
									_struct1.getText(), _struct2.getText());
							_vp.showRNAInterpolated(getRNA1(), m);
						} catch (MappingException e3) {
							_vp.drawRNAInterpolated(_seq.getText(),
									getStruct1());
						}
					} else {
						try {
							Mapping m = Mapping.readMappingFromAlignment(
									_struct2.getText(), _struct1.getText());
							_vp.showRNAInterpolated(getRNA2(), m);
						} catch (MappingException e3) {
							_vp.drawRNAInterpolated(_seq.getText(),
									getStruct2());
						}
					}
					_vp.repaint();
				} catch (ExceptionNonEqualLength e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}
			}
		});

		_seqPanel.setLayout(new BorderLayout());
		_seqPanel.add(_seqLabel, BorderLayout.WEST);
		_seqPanel.add(_seq, BorderLayout.CENTER);

		_struct1Label.setPreferredSize(new Dimension(marginTools, 15));
		_struct1Label.setHorizontalTextPosition(JLabel.LEFT);
		_struct1.setFont(textFieldsFont);
		_struct1Panel.setLayout(new BorderLayout());
		_struct1Panel.add(_struct1Label, BorderLayout.WEST);
		_struct1Panel.add(_struct1, BorderLayout.CENTER);

		_struct2Label.setPreferredSize(new Dimension(marginTools, 15));
		_struct2Label.setHorizontalTextPosition(JLabel.LEFT);
		_struct2.setFont(textFieldsFont);
		_struct2Panel.setLayout(new BorderLayout());
		_struct2Panel.add(_struct2Label, BorderLayout.WEST);
		_struct2Panel.add(_struct2, BorderLayout.CENTER);

		_input.setLayout(new GridLayout(3, 0));
		_input.add(_seqPanel);
		_input.add(_struct1Panel);
		_input.add(_struct2Panel);

		JPanel goPanel = new JPanel();
		goPanel.setLayout(new BorderLayout());

		_tools.setLayout(new BorderLayout());
		_tools.add(_input, BorderLayout.CENTER);
		_tools.add(_info, BorderLayout.SOUTH);
		_tools.add(goPanel, BorderLayout.EAST);

		goPanel.add(_goButton, BorderLayout.CENTER);
		goPanel.add(_switchButton, BorderLayout.SOUTH);

		getContentPane().setLayout(new BorderLayout());
		getContentPane().add(_vp, BorderLayout.CENTER);
		getContentPane().add(_tools, BorderLayout.SOUTH);

		setVisible(true);
		_vp.getVARNAUI().UIRadiate();
	}

	public RNA getRNA1() {

		if (!_str1Backup.equals(getStruct1())) {
			try {
				_RNA1.setRNA(_seq.getText(), getStruct1());
				_RNA1.drawRNA(_vp.getDrawMode());
			} catch (ExceptionUnmatchedClosingParentheses e) {
				e.printStackTrace();
			} catch (ExceptionFileFormatOrSyntax e1) {
				_vp.errorDialog(e1);
			} catch (ExceptionDrawingAlgorithm e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			_str1Backup = getStruct1();
		}
		return _RNA1;
	}

	public RNA getRNA2() {
		if (!_str2Backup.equals(getStruct2())) {
			try {
				_RNA2.setRNA(_seq.getText(), getStruct2());
				_RNA2.drawRNA(_vp.getDrawMode());
			} catch (ExceptionUnmatchedClosingParentheses e) {
				e.printStackTrace();
			} catch (ExceptionFileFormatOrSyntax e1) {
				_vp.errorDialog(e1);
			} catch (ExceptionDrawingAlgorithm e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			_str2Backup = getStruct2();
		}
		return _RNA2;
	}

	public String getStruct1() {
		return cleanStruct(_struct1.getText());
	}

	public String getStruct2() {
		return cleanStruct(_struct2.getText());
	}

	private String cleanStruct(String struct) {
		struct = struct.replaceAll("[:-]", "");
		return struct;
	}

	public String[][] getParameterInfo() {
		String[][] info = {
				// Parameter Name Kind of Value Description,
				{ "sequenceDBN", "String", "A raw RNA sequence" },
				{ "structureDBN", "String",
						"An RNA structure in dot bracket notation (DBN)" },
				{ errorOpt, "boolean", "To show errors" }, };
		return info;
	}

	public void init() {
		_vp.setBackground(_backgroundColor);
		_error = true;
	}

	@SuppressWarnings("unused")
	private Color getSafeColor(String col, Color def) {
		Color result;
		try {
			result = Color.decode(col);
		} catch (Exception e) {
			try {
				result = Color.getColor(col, def);
			} catch (Exception e2) {
				return def;
			}
		}
		return result;
	}

	public VARNAPanel get_varnaPanel() {
		return _vp;
	}

	public void set_varnaPanel(VARNAPanel surface) {
		_vp = surface;
	}

	public JTextField get_struct() {
		return _struct1;
	}

	public void set_struct(JTextField _struct) {
		this._struct1 = _struct;
	}

	public JTextField get_seq() {
		return _seq;
	}

	public void set_seq(JTextField _seq) {
		this._seq = _seq;
	}

	public JLabel get_info() {
		return _info;
	}

	public void set_info(JLabel _info) {
		this._info = _info;
	}

	public static void main(String[] args) {
		SuperpositionDemo d = new SuperpositionDemo();
		d.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		d.pack();
		d.setVisible(true);
	}
}
