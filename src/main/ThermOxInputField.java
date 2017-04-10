package main;

import javax.swing.*;
import java.text.DecimalFormat;

public class ThermOxInputField extends JTextField {
    static final DecimalFormat DECIMAL_FORMATTER_DEFAULT = new DecimalFormat("0.0000");

    private String labelName;   // optional
    private JLabel jLabel;  // optional
    private double defaultPresetValue;  // required
    private boolean isUserEditable = true;  // optional

    public static class Builder {
        private String labelName;
        private JLabel jLabel;
        private double defaultPresetValue;
        private boolean isUserEditable = true;

        public Builder(double defaultPresetValue) {
            this.defaultPresetValue = defaultPresetValue;
        }

        public Builder labelName (String labelName) {
            this.labelName = labelName;
            return this;
        }

        public Builder jLabel (JLabel jLabel) {
            this.jLabel = jLabel;
            return this;
        }

        public Builder userEditable (boolean isUserEditable) {
            this.isUserEditable = isUserEditable;
            return this;
        }

        public ThermOxInputField build() {
            return new ThermOxInputField(this);
        }
    }

    private ThermOxInputField(Builder builder) {
        this.labelName = builder.labelName;
        this.jLabel = builder.jLabel;
        this.defaultPresetValue = builder.defaultPresetValue;
        this.isUserEditable = builder.isUserEditable;

        this.setEditable(isUserEditable);
        this.setFocusable(isUserEditable);

        if (!isUserEditable) {
            this.setBackground(new java.awt.Color(230, 230, 230));
        }

        this.setText("0.0");
        this.setHorizontalAlignment(javax.swing.JTextField.TRAILING);

        if (this.jLabel != null && this.labelName != null) {
            this.jLabel.setText(this.labelName);
        }
    }

    public double getDefaultPresetValue() {
        return this.defaultPresetValue;
    }

    public String getLabelName() {
        return this.labelName;
    }

    public JLabel getjLabel() {
        return this.jLabel;
    }

    public boolean isUserEditable() {
        return this.isUserEditable;
    }

    public void resetValue () {
        this.setText(DECIMAL_FORMATTER_DEFAULT.format(this.defaultPresetValue));
    }

    public void resetValue (boolean useDefaultFormatter) {
        String text = useDefaultFormatter ? DECIMAL_FORMATTER_DEFAULT.format(this.defaultPresetValue) : String.valueOf(this.defaultPresetValue);
        this.setText(text);
    }
}