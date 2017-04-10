package main;

import javax.swing.*;
import java.text.DecimalFormat;

public class ThermalOxidizerInputTextField extends JTextField {
    static final DecimalFormat DECIMAL_FORMATTER_DEFAULT = new DecimalFormat("0.0000");

    private double defaultValue = 0.0;
    private double clearVal = 0.0;
    private String fieldName;
    private JLabel jLabel;

    private boolean isUserEditable = true;

    public ThermalOxidizerInputTextField (double defaultValue, boolean isUserEditable) {
        this(null, null, defaultValue, isUserEditable);
    }

    public ThermalOxidizerInputTextField (String fieldName, double defaultValue) {
        this(fieldName, null, defaultValue, true);
    }

    public ThermalOxidizerInputTextField (String fieldName, JLabel jLabel, double defaultValue) {
        this(fieldName, jLabel, defaultValue, true);
    }

    public ThermalOxidizerInputTextField (String fieldName, double defaultValue, boolean isUserEditable) {
        this(fieldName, null, defaultValue, isUserEditable);
    }

    public ThermalOxidizerInputTextField (String fieldName, JLabel jLabel, double defaultValue, boolean isUserEditable) {
        this.fieldName = fieldName;
        this.defaultValue = defaultValue;
        this.isUserEditable = isUserEditable;

        this.setEditable(isUserEditable);
        this.setFocusable(isUserEditable);

        if (!isUserEditable) {
            this.setBackground(new java.awt.Color(230, 230, 230));
        }

        this.setText("0.0");
        this.setHorizontalAlignment(javax.swing.JTextField.TRAILING);

        if (jLabel != null) {
            this.jLabel = jLabel;
            this.jLabel.setText(fieldName);
        }
    }

    public void clearValue () {
        this.setText(DECIMAL_FORMATTER_DEFAULT.format(this.clearVal));
    }

    public double getDefaultValue() {
        return this.defaultValue;
    }

    public String getFieldName() {
        return this.fieldName;
    }

    public JLabel getjLabel() {
        return this.jLabel;
    }

    public boolean isUserEditable() {
        return this.isUserEditable;
    }

    public void resetValue () {
        this.setText(DECIMAL_FORMATTER_DEFAULT.format(this.defaultValue));
    }

    public void resetValue (boolean useDefaultFormatter) {
        String text = useDefaultFormatter ? DECIMAL_FORMATTER_DEFAULT.format(this.defaultValue) : String.valueOf(this.defaultValue);
        this.setText(text);
    }

    public void setDefaultValue(double defaultValue) {
        this.defaultValue = defaultValue;
    }

    public void setFieldName(String fieldName) {
        this.fieldName = fieldName;
    }

    public void setjLabel(JLabel jLabel) {
        this.jLabel = jLabel;
    }

    public void setUserEditable(boolean userEditable) {
        isUserEditable = userEditable;
    }
}
